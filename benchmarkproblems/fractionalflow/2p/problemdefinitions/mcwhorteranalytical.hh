#ifndef DUNE_BUCKLEYLEVERETTANALYTICAL_HH
#define DUNE_BUCKLEYLEVERETTANALYTICAL_HH

#include"dumux/fractionalflow/exsolution.hh"

/**
 * @file
 * @brief Class adding the analytical solution of the Buckley Leverett problem
 * @author Markus Wolff
 */

namespace Dune
{

	template<class G, class RT, class VC> class McWWithAnalytical : public ExSolution<G, RT>,
		public McWhorterTransportProblem<G, RT,VC>
	{
		enum
		{	n=G::dimension,dimworld = G::dimensionworld,m=2};
		typedef typename G::ctype DT;

		typedef typename G::LevelGridView GV;
	    typedef typename GV::IndexSet IS;
		typedef typename GV::template Codim<0>::Iterator Iterator;

		typedef BlockVector<FieldVector<RT, m> > BVu;
		typedef BlockVector<FieldVector<RT, 1> > BV;

public:

		McWWithAnalytical(VC& variables, TwoPhaseRelations &law = *(new LinearLaw), const FieldVector<DT,n> Left = 0,
				const FieldVector<DT,n> Right = 0, RT cfl = 1, RT permeability = 1e-10) :
			ExSolution<G, RT>(variables.grid), McWhorterTransportProblem<G, RT, VC>(
					variables, law, Left, Right, true), CFL(cfl), K(permeability),
					tolanalytic(1e-14), R(0)
		{
			prepareanalytical();
		}

		void prepareanalytical()
		{
			Swr=this->materialLaw.wettingPhase.Sr();
			Snr=this->materialLaw.nonwettingPhase.Sr();
			Porosity=this->porosity();
			Sinit=this->Sinit_;

			time=0;

			h= (1-Sinit)/intervalnum;
			//std::cout<<"h = "<<h_<<std::endl;

			//define saturation range for analytic solution
			SatVec= 0;
			for (int i=1; i<pointnum; i++)
			{
				SatVec[i]=SatVec[i-1]+h;
			}

			//get fractional flow function vector
			for (int i=0; i<pointnum; i++)
			{
				fractionalW[i] = this->materialLaw.fractionalW(SatVec[i]);
			}

			//get capillary pressure derivatives
			dpcdsw=0;
			for (int i=0; i<pointnum; i++)
			{
				dpcdsw[i] = this->materialLaw.dPdS(SatVec[i]);
			}

			//set initial fW
			if (Sinit == 0)
				finit=0;
			else
				finit=fractionalW[0];

			fractionalW[0]=0;

			//normalize fW
			for (int i=0; i<pointnum; i++)
			{
				fn[i]= R * (fractionalW[i] - finit)/ (1 - R * finit);
			}

			//diffusivity function
			for (int i=0; i<pointnum; i++)
			{
				D[i] = fractionalW[i]*this->materialLaw.mobN(1-SatVec[i])*(-dpcdsw[i])*K;
			}

			for (int i=0; i<pointnum; i++)
			{
				Gk[i] = D[i]/(1-fn[i]);
			}
			Gk[0] = 0;
			//std::cout<<"Gk_ = "<<Gk_<<std::endl;
			return;
		}

		void settime(double &dt)
		{
			time+=dt*CFL;
			return;
		}

		void updateExSol()
		{
			double Ak = 0;
			double Akm1 = 0;
			double diff = 1e100;
			int k = 0;
			a=0;
			b=0;
			Fp=0;
			double I0 = 0;
			double Ii = 0;

			while (diff > tolanalytic)
			{
				k++;
				//std::cout<<" k = "<<k<<std::endl;
				if (k > 50000)
				{
					std::cout<<"Analytical solution: Too many iteratioins!"<<std::endl;
					break;
				}

				Akm1=Ak;
				I0=0;
				for (int i=0; i<intervalnum; i++)
				{
					a[i] = 0.5 * h * Sinit *(Gk[i] + Gk[i+1])+ pow(h, 2) / 6* ((3* i + 1) * Gk[i]
							+ (3 * i + 2) * Gk[i+1]);
					b[i] = 0.5 * h * (Gk[i] + Gk[i+1]);
					I0 += (a[i] - Sinit * b[i]);
				}
				//std::cout<<" I0 = "<<I0<<std::endl;

				Gk[0]=0;
				for (int i=1; i<pointnum; i++)
				{
					Ii=0;
					for (int j = i; j<intervalnum; j++)
						Ii += (a[j] - SatVec[i] * b[j]);
					//Gk[i] = D[i] + Gk[i]*(fn[i] + Ii/I0); // method A  
					Gk[i] = (D[i] + Gk[i]*fn[i])/(1 - Ii/I0); // method B
				}

				Ak = pow((0.5*Porosity/pow((1 - finit), 2)*I0), 0.5);
				diff=fabs(Ak - Akm1);
				//std::cout<<"diff = "<<diff<<std::endl;
			}
			//  std::cout<<" b = "<<b<<std::endl;
			//  std::cout<<" Ak = "<<Ak<<std::endl;
			for (int i = 0; i<pointnum ; i++)
			{
				for (int j = i; j<intervalnum; j++)
					Fp[i] += b[j]/I0;
			}
			//std::cout<<" Fp = "<<Fp<<std::endl;
			for (int i = 0; i<pointnum ; i++)
			{
				xf[i]= 2 * Ak * (1 - finit * R)/ Porosity * Fp[i]* pow(time, 0.5);
			}
			//std::cout<<" xf = "<<xf_<<std::endl;  

			//iterate over vertices and get analytical saturation solution

			const GV& gridview(this->grid.levelView(0));

			Iterator eendit = gridview.template end<0>();
			for (Iterator it = gridview.template begin<0>(); it!= eendit; ++it)
			{

				//find index of current vertex
				int index = this->mapper.map(*it);

				//get global coordinate of vertex
				Dune::FieldVector<DT,dimworld> global = it->geometry()[0];

				//find x_f next to global coordinate of the vertex
				int xnext = 0;
				for (int i=intervalnum; i>0; i--)
				{
					if (global[0]<xf[i])
					{
						xnext = i;
						break;
					}
				}
				//account for the area not yet reached by the front

				if (global[0] > xf[0])
				{
					this->uExInVertex(Sinit, index, 0);
					continue;
				}

				this->uExInVertex(SatVec[xnext], index, 0);
				//std::cout<<"Analytical = "<<SatVec_[xnext]<<std::endl;
			}

			//call function to calculate the saturation error
			this->calcSatError(this->variables.saturation);

			return;
		}

		BlockVector<FieldVector<RT, 2> >& getuEx()
		{
			return this->uEx;
		}

private:
		RT Swr;
		RT Snr;
		RT Porosity;
		RT Sinit;

		RT time;
		RT CFL;
		RT K;
		RT tolanalytic;
		enum
		{	intervalnum = 1000, pointnum = intervalnum+1};
		FieldVector<RT, pointnum> SatVec;
		FieldVector<RT,pointnum> fractionalW;
		FieldVector<RT, pointnum> dpcdsw;
		FieldVector<RT, pointnum> fn;
		FieldVector<RT, pointnum> D;
		FieldVector<RT, pointnum> Gk;
		FieldVector<RT, pointnum> xf;
		RT finit;
		RT R;
		RT h;
		FieldVector<RT,intervalnum> a;
		FieldVector<RT,intervalnum> b;
		FieldVector<RT,pointnum> Fp;
	};

}
#endif
