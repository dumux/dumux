#ifndef DUNE_CO2HEATFLUX_HH
#define DUNE_CO2HEATFLUX_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/linearlaw.hh>
#include<dumux/material/brookscoreylaw.hh>
#include<dumux/2penergy/2penergyproblem.hh>


/**
 * @file
 * @brief  Base class for defining an instance of the TwoPhase problem
 * @author Bernd Flemisch
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
  class CO2HeatFlux : public TwoPhaseHeatProblem<G, RT> {
	typedef typename G::ctype DT;
	enum {n=G::dimension, m=3};
	enum {swrIdx=0, snrIdx=1, lamIdx=2, pbIdx=3};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
	virtual const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi)
	{
		return permloc_;
	}

	virtual FieldVector<RT,m> q (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	{
		FieldVector<RT,m> values(0);

		return values;
	}

	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,n>& x, const Entity& e, 
					const IntersectionIterator& intersectionIt, 
					   const FieldVector<DT,n>& xi) const 
	{
		FieldVector<BoundaryConditions::Flags, m> values(Dune::BoundaryConditions::dirichlet); 
		if(x[0]<1.e-2 && x[1] > 1. && x[1] < 3.)
		{
			values[0] = Dune::BoundaryConditions::neumann;
			values[1] = Dune::BoundaryConditions::neumann;
		        values[2] = Dune::BoundaryConditions::neumann;
		}
		return values;
	}
	
	virtual void dirichletIndex(const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi, FieldVector<int,m>& dirichletIndex) const 
	{
		for (int i = 0; i < m; i++)
			dirichletIndex[i]=i;
		return;
	}

	virtual FieldVector<RT,m> g (const FieldVector<DT,n>& x, const Entity& e, 
				const IntersectionIterator& intersectionIt, 
				  const FieldVector<DT,n>& xi) const 
	{
		FieldVector<RT,m> values(0);
		
		values[0] = p0_; 
		values[1] = 0.3;
		values[2] = 311.; 
		
		return values;
	}
	  
	virtual FieldVector<RT,m> J (const FieldVector<DT,n>& x, const Entity& e, 
				const IntersectionIterator& intersectionIt, const FieldVector<DT,n>& xi) const 
	{
		FieldVector<RT,m> values(0);

		values[0] = 0.;
		values[1] = 0.;
		values[2] = 10000.77;
		
		return values;
	}
	
	// Initial Conditions for global vector x, element e and local vector xi
	virtual FieldVector<RT,m> initial (const FieldVector<DT,n>& x, const Entity& e, 
				  const FieldVector<DT,n>& xi) const 
	{
		FieldVector<RT,m> values(0);

		values[0] = p0_; 
		values[1] = 0.3;
		values[2] = 311.;
	
		return values;
	}

	double porosity (const FieldVector<DT,n>& x, const Entity& e, 
			  const FieldVector<DT,n>& xi) const 
	{
		return 0.2;
	}
	
	virtual FieldVector<RT,4> soilParameters (const FieldVector<DT,n>& x, const Entity& e, 
			  const FieldVector<DT,n>& xi) const
	{
		FieldVector<RT,4> values(0);
		enum{soilDensity = 0,soilHeatCap = 1,soilLambdaDry = 2, soilLambdaSw = 3};
		
		values[soilDensity] = 2650.0;
		values[soilHeatCap] = 800.0;
		values[soilLambdaDry] = 0.32;
		values[soilLambdaSw] = 2.7;
		
		return values;
	}
	
	virtual FieldVector<RT,n> gravity () const 
	{
		FieldVector<RT,n> values(0);
		
		values[1] = -9.81;
		
		return values;
	}
	  
	virtual FieldVector<RT,4> materialLawParameters (const FieldVector<DT,n>& x, const Entity& e, 
			  const FieldVector<DT,n>& xi) const 
	{
		FieldVector<RT,4> values(0);
		values[swrIdx] = swr_;
		values[snrIdx] = snr_;
		values[lamIdx] = lambda_;
		values[pbIdx] = pb_;
		
		return values;
	}

	CO2HeatFlux(TwoPhaseRelations& law = *(new BrooksCoreyLaw), RT pdown = 3.086e5, RT swr = 0.2, RT snr = 0.05, RT pb = 10000.0, RT lambda = 2.0 ) 
	: TwoPhaseHeatProblem<G, RT>(law) 
	{	
		p0_ = pdown;
		swr_ = swr;
		snr_ = snr;
		pb_ = pb;
		lambda_ = lambda;
		permloc_ = 0; 
		
		for (int i = 0; i < n; i++)
			permloc_[i][i] = 1.0e-14;


	}
	
	private:
		Dune::FieldMatrix<DT,n,n> permloc_;
		Dune::FieldMatrix<DT,n,n> permlocWell_;
		Dune::FieldMatrix<DT,n,n> permlocAquitard_;
		RT p0_;
		RT swr_, snr_;
		RT pb_, lambda_;
		RT soilDens_, soilHeatCp_, soilLDry_, soilLSw_;
		
  };

}
#endif
