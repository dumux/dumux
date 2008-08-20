#ifndef DUNE_WATERCO2PROBLEM_HH
#define DUNE_WATERCO2PROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations_deprecated.hh>
#include<dumux/2p2cni/2p2cniproblem.hh>
#include<dumux/material/linearlaw_deprecated.hh>
#include<dumux/material/brookscoreylaw_deprecated.hh>


/**
 * @file
 * @brief  Base class for defining an instance of the TwoPhase problem
 * @author Bernd Flemisch, Melanie Darcis, Klaus Mosthaf
 */

namespace Dune
{
  //! base class that defines the parameters of a diffusion equation
  /*! An interface for defining parameters for the stationary diffusion equation 
   * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$, 
   * \f$p = g\f$ on \f$\Gamma_1\f$, and \f$\lambda K \text{grad}\, p = J\f$ 
   * on \f$\Gamma_2\f$. Here, 
   * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability, 
   * and \f$\lambda\f$ the total mobility, possibly depending on the 
   * saturation. 
   *
   *	Template parameters are:
   *	
   *	- Grid  a DUNE grid type
   *	- RT    type used for return values 
   */
  template<class G, class RT>
  class WaterCO2Problem : public TwoPTwoCNIProblem<G, RT> {
	typedef typename G::ctype DT;
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;
	enum {dim=G::dimension, m=3};
	enum {swrIdx=0, snrIdx=1, lamIdx=2, pbIdx=3};
	enum {wPhase = 0, nPhase = 1};
	enum {pWIdx = 0, switchIdx = 1, teIdx = 2};


  public:

	// permeabilities
	virtual const FieldMatrix<DT,dim,dim>& K (const FieldVector<DT,dim>& x, const Entity& e, 
					const FieldVector<DT,dim>& xi)
	{
		return permloc_;
	}

	virtual const FieldMatrix<DT,dim,dim>& K (const FieldVector<DT,dim>& x)
	//, const Entity& e, const FieldVector<DT,dim>& xi)
	{
		return permloc_;
	}
		
	virtual const FieldVector<RT,2>& D (const FieldVector<DT,dim>& x) const 
	{
		return diffusion_;
	}
		
	virtual FieldVector<RT,m> q (const FieldVector<DT,dim>& x, const Entity& e, 
					const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		return values;
	}

/////////////////////////////
// TYPE of the boundaries
/////////////////////////////
	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,dim>& x, const Entity& e, 
					const IntersectionIterator& intersectionIt, 
					   const FieldVector<DT,dim>& xi) const 
	{
		// initialize all boundaries as Neumann
		FieldVector<BoundaryConditions::Flags, m> values(Dune::BoundaryConditions::neumann); 

		if(x[0] < 1e-10)
		{
			values = Dune::BoundaryConditions::dirichlet;
		}
		
//		if(x[1] < 1e-10)
//		{
//			values = Dune::BoundaryConditions::dirichlet;
//		}

		return values;
	}
	
	virtual void dirichletIndex(const FieldVector<DT,dim>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,dim>& xi, FieldVector<int,m>& dirichletIndex) const 
	{
		for (int i = 0; i < m; i++)
			dirichletIndex[i]=i;
		return;
	}

/////////////////////////////
// DIRICHLET boundaries
/////////////////////////////
	virtual FieldVector<RT,m> g (const FieldVector<DT,dim>& x, const Entity& e, 
				const IntersectionIterator& intersectionIt, 
				  const FieldVector<DT,dim>& xi) const 
	{
		FieldVector<RT,m> values(0.);
		
		values[pWIdx] = 1e5 + (depthBOR_ - x[1])*1000*9.81; 
		values[switchIdx] = 0.;
		values[teIdx] = 283. + (depthBOR_ - x[1])*0.03;

//		RT pinj = 0.0;	
//		if (x[0] < 1e-10 && x[1] > 1.0 && x[1] < 2.0)
//		{
//			values[pWIdx] = 1e5 + (depthBOR_ - x[1])*1000*9.81 + pinj; 
//			values[switchIdx] = 0.3;
//			values[teIdx] = 283. + (depthBOR_ - x[1])*0.03;	
//		}
			
		return values;
	}
	
/////////////////////////////
// NEUMANN boundaries
/////////////////////////////	  
	virtual FieldVector<RT,m> J (const FieldVector<DT,dim>& x, const Entity& e, 
				const IntersectionIterator& intersectionIt, const FieldVector<DT,dim>& xi) const 
	{
		FieldVector<RT,m> values(0.0);
		
		// negative values for injection		
		if (x[1] > 1.0 && x[1] < 5.0)
		{
			values[pWIdx] = 0.0;
			values[switchIdx] = -1e-5;
			values[teIdx] = 0.0;
		}
		
		return values;
	}
	
/////////////////////////////
// INITIAL values
/////////////////////////////
	virtual FieldVector<RT,m> initial (const FieldVector<DT,dim>& x, const Entity& e, 
				  const FieldVector<DT,dim>& xi) const 
	{
		FieldVector<RT,m> values(0);

		values[pWIdx] = 1e5 + (depthBOR_ - x[1])*1000*9.81; 
		values[switchIdx] = 0.;
		values[teIdx] = 283. + (depthBOR_ - x[1])*0.03;
	
		return values;
	}

	int initialPhaseState (const FieldVector<DT,dim>& x, const Entity& e, 
				  const FieldVector<DT,dim>& xi) const 
	{

		enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase states
		int state;

//		if (x[1] >= innerLowerLeft_[1] && x[1] <= innerUpperRight_[1]
//		      && x[0] >= innerLowerLeft_[0])
		state = waterPhase;
			
		return state;
	}
//////////////////////////////	  

	double porosity (const FieldVector<DT,dim>& x, const Entity& e, 
			  const FieldVector<DT,dim>& xi) const 
	{
		return 0.2;
	}
	
	virtual FieldVector<RT,4> soilParameters (const FieldVector<DT,dim>& x, const Entity& e, 
			  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,4> values(0);
		enum{soilDensity = 0,soilHeatCap = 1,soilLambdaDry = 2, soilLambdaSw = 3};
		
		values[soilDensity] = 2650.0;
		values[soilHeatCap] = 800.0;
		values[soilLambdaDry] = 0.32;
		values[soilLambdaSw] = 2.7;
		
		return values;
	}
	
	virtual FieldVector<RT,dim> gravity () const 
	{
		FieldVector<RT,dim> values(0);
		
		values[1] = -9.81;
		
		return values;
	}
	  
	virtual FieldVector<RT,4> materialLawParameters (const FieldVector<DT,dim>& x, const Entity& e, 
			  const FieldVector<DT,dim>& xi) const 
	{
		FieldVector<RT,4> values(0);
		values[swrIdx] = swr_;
		values[snrIdx] = snr_;
		values[lamIdx] = lambda_;
		values[pbIdx] = pb_;
		
		return values;
	}

	WaterCO2Problem(TwoPhaseRelations& law = *(new BrooksCoreyLaw), 
			MultiComp& multicomp = *(new CWaterAir), RT depthBOR = 0.0, 
			RT swr = 0.1, RT snr = 0.1, RT pb = 10000.0, RT lambda = 2.0) 
	: TwoPTwoCNIProblem<G, RT>(law, multicomp) 
	{	
		depthBOR_ = depthBOR;
		swr_ = swr;
		snr_ = snr;
		pb_ = pb;
		lambda_ = lambda;
		permloc_ = 0; 
		
		for (int i = 0; i < dim; i++)
			permloc_[i][i] = 9.2e-12;
		
		diffusion_[wPhase] = 2.0E-5; // diffusion coefficient for water in gas phase
		diffusion_[nPhase] = 2.6E-9; // diffusion coefficient for co2 in water phase

	}
	
	private:
		Dune::FieldMatrix<DT,dim,dim> permloc_;
		Dune::FieldVector<RT,2> diffusion_;
		RT depthBOR_;
		RT swr_, snr_;
		RT pb_, lambda_;
		RT soilDens_, soilHeatCp_, soilLDry_, soilLSw_;
		
  };

}
#endif
