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

    template<class G, class RT, class VC> class BLWithAnalytical :
        public ExSolution<G, RT>,
        public BuckleyLeverettTransportProblem<G, RT,VC>
    {
        enum
        {n=G::dimension,dimworld = G::dimensionworld,m=2};
        typedef typename G::ctype DT;
//        typedef typename G::template Codim <0>:: LeafIterator Iterator;

        typedef typename G::LevelGridView GV;
        typedef typename GV::IndexSet IS;
        typedef typename GV::template Codim<0>::Iterator Iterator;

        typedef BlockVector<FieldVector<RT, m> > BVu;
        typedef BlockVector<FieldVector<RT, 1> > BV;

public:

        BLWithAnalytical(VC& variables, TwoPhaseRelations &law = *(new LinearLaw), const FieldVector<DT,n> Left = 0,
                const FieldVector<DT,n> Right = 0, RT cfl = 1,
                RT totalvelocity = 3e-7) :
            ExSolution<G, RT>(variables.grid), BuckleyLeverettTransportProblem<
                    G, RT, VC>(variables, law, Left, Right, true), CFL(cfl),
                    vtot(totalvelocity)
        {
            prepareanalytical();
        }

        void prepareanalytical()
        {
            Swr = this->materialLaw.wettingPhase.Sr();
            Snr = this->materialLaw.nonwettingPhase.Sr();

            time=0;

            SatVec=Swr;
            for (int i=1; i<pointnum; i++)
            {
                SatVec[i]=SatVec[i-1]+(1-Snr-Swr)/intervalnum;
            }
            //        std::cout<<"SatVec = "<<SatVec<<std::endl;
            for (int i=0; i<pointnum; i++)
            {
                fractionalW[i] = this->materialLaw.fractionalW(SatVec[i]);
            }
            //        std::cout<<"fractionalW = "<<fractionalW<<std::endl;
            dfwdsw=0;
            for (int i=1; i<intervalnum; i++)
            {
                dfwdsw[i]=(fractionalW[i+1]-fractionalW[i-1])/(SatVec[i+1]
                        -SatVec[i-1]);
            }
            //        std::cout<<"dfwdsw = "<<dfwdsw<<std::endl;
            for (int i=0; i<pointnum; i++)
            {
                if (dfwdsw[i]>dfwdsw[i+1])
                {
                    dfwdswmax = i;
                    break;
                }
            }
            return;
        }

        void settime(double &dt)
        {
            time+=dt*CFL;
            return;
        }

        void updateExSol()
        {
            //possition of the fluid front
            xf=0;
            for (int i=0; i<pointnum; i++)
            {
                xf[i]=vtot*time/this->porosity()*dfwdsw[i];
            }
            //std::cout<<"xf = "<<xf_<<std::endl;
            int xhelp=pointnum/3;
            int xhelpold = 0;
            int xhelpoldold = 0;
            int xfmax = 0;
            ;

            //position of maximum xf
            for (int i=0; i<pointnum; i++)
            {
                if (xf[i]>xf[i+1])
                {
                    xfmax=i;
                    break;
                }
            }

            //balancing of the areas ahead of the front and below the curve
            bool a = true;
            double A1;
            double A2;
            double b;
            int xhelp2=0;

            while (a)
            {
                if (this->materialLaw.isLinear())
                    break;

                A1=0;
                for (int i=0; i<=xhelp-1; i++)
                {
                    A1+=(SatVec[i]-Swr+SatVec[i+1]-Swr)*0.5*(xf[i+1]-xf[i]);
                }
                A2=0;
                for (int i=xhelp; i<=xfmax-1; i++)
                {
                    A2+=(SatVec[xfmax]-SatVec[i]+SatVec[xfmax]-SatVec[i+1])*0.5
                            *(xf[i+1]-xf[i]);
                }
                b=xf[xfmax];
                xhelp2=xfmax;
                while (b>xf[xhelp])
                {
                    xhelp2+=1;
                    b=xf[xhelp2];
                }
                for (int i=xfmax; i<=xhelp2; i++)
                {
                    A2+=(SatVec[i]-SatVec[xfmax]+SatVec[i+1]-SatVec[xfmax])*0.5
                            *(xf[i]-xf[i+1]);
                }
                xhelpoldold=xhelpold;
                xhelpold=xhelp;
                if (fabs(A1)>fabs(A2))
                {
                    xhelp = xhelp - 1;
                }
                if (fabs(A1)<fabs(A2))
                {
                    xhelp = xhelp +1;
                }
                if (xhelp == xhelpoldold)
                {
                    a=false;
                }
            }
            //        std::cout<<"xf = "<<xf<<std::endl;
            //iterate over vertices and get analytical saturation solution

            const GV& gridview(this->grid.levelView(0));

            Iterator eendit = gridview.template end<0>();
            for (Iterator it = gridview.template begin<0>(); it
                    != eendit; ++it)
            {

                //find index of current vertex
                int index = this->mapper.map(*it);

                //get global coordinate of vertex
                Dune::FieldVector<DT,dimworld> global = it->geometry()[0];

                //account for linear material law
                if (this->materialLaw.isLinear())
                {
                    if (global[0]<=xf[1])
                    {
                        RT oneminusSnr = 1-Snr;
                        this->uExInVertex(oneminusSnr, index, 0);
                    }
                    if (global[0]>xf[1])
                        this->uExInVertex(Swr, index, 0);
                }

                //non-linear material law
                else
                {
                    //find x_f next to global coordinate of the vertex
                    int xnext = 0;
                    for (int i=intervalnum; i>=0; i--)
                    {
                        if (global[0]<xf[i])
                        {
                            xnext = i;
                            break;
                        }
                    }

                    //account for the area not yet reached by the front
                    if (global[0] > xf[xhelp2])
                    {
                        this->uExInVertex(Swr, index, 0);
                        continue;
                    }

                    if (global[0] <= xf[xhelp2])
                    {
                        this->uExInVertex(SatVec[xnext], index, 0);
                        continue;
                    }
                }
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

        RT time;
        RT CFL;
        RT vtot;
        enum{intervalnum = 1000, pointnum = intervalnum+1};
        FieldVector<RT, pointnum> SatVec;
        FieldVector<RT,pointnum> fractionalW;
        FieldVector<RT,pointnum > dfwdsw;
        FieldVector<RT,pointnum> xf;
        int dfwdswmax;
    };

}
#endif
