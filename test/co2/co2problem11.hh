#ifndef DUNE_CO2PROBLEM11_HH
#define DUNE_CO2PROBLEM11_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations_deprecated.hh>
#include<dumux/material/linearlaw_deprecated.hh>
#include<dumux/twophase/twophaseproblem.hh>

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
  class CO2Problem11 : public TwoPhaseProblem<G, RT> {
	typedef typename G::ctype DT;
	enum {n=G::dimension, m=2};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
	virtual const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi)
	{
		if (x[0] > -0.132934 && x[0] < 0.132934
				&& x[1] > -0.132934 && x[1] < 0.132934)
			return permlocWell;
	
		return permloc;
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
		FieldVector<BoundaryConditions::Flags, m> values(Dune::BoundaryConditions::neumann); 

	        if (x[0] < -500+1e-3 || x[0] > 500-1e-3 || x[1] < -500+1e-3 || x[1] > 500-1e-3) 
	        	values = Dune::BoundaryConditions::dirichlet;
		
		return values;
	}

	virtual FieldVector<RT,m> g (const FieldVector<DT,n>& x, const Entity& e, 
				const IntersectionIterator& intersectionIt, 
				  const FieldVector<DT,n>& xi) const 
	{
		FieldVector<RT,m> values(0);
		values[0] = p0 - 1045.0*9.81*x[2]; 
		values[1] = 0.0;

		return values;
	}
	  
	virtual FieldVector<RT,m> J (const FieldVector<DT,n>& x, const Entity& e, 
				const IntersectionIterator& intersectionIt, 
				  const FieldVector<DT,n>& xi) const 
	{
		FieldVector<RT,m> values(0);
		if ((x[0]+100.0)*(x[0]+100.0) + x[1]*x[1] < 0.3 && x[2] < 30.0-1e-3 && x[2] > 1e-3) 
			values[1] = -0.27802;
		
		return values;
	}
	  
	virtual FieldVector<RT,m> initial (const FieldVector<DT,n>& x, const Entity& e, 
				  const FieldVector<DT,n>& xi) const 
	{
		FieldVector<RT,m> values(0);
		values[0] = p0 - 1045.0*9.81*x[2];
		values[1] = 0.0;
	
		return values;
	}

	double porosity (const FieldVector<DT,n>& x, const Entity& e, 
			  const FieldVector<DT,n>& xi) const 
	{
		return 0.15;
	}
	
	virtual FieldVector<RT,n> gravity () const 
	{
		FieldVector<RT,n> values(0);
		
		values[2] = -9.81;
		
		return values;
	}
	  
	virtual FieldVector<RT,4> materialLawParameters (const FieldVector<DT,n>& x, const Entity& e, 
			  const FieldVector<DT,n>& xi) const 
	{
		FieldVector<RT,4> values(0);
		
		return values;
	}

	CO2Problem11(TwoPhaseRelations& law = *(new LinearLaw), RT pdown = 3.086e7) 
	: TwoPhaseProblem<G, RT>(law) 
	{	
		p0 = pdown;
		permloc = 0; 
		permlocWell = 0;
		
		for (int i = 0; i < n; i++)
			permloc[i][i] = 2.0e-14;
		for (int i = 0; i < n; i++)
			permlocWell[i][i] = 1.0e-12;


	}
	
	private:
		Dune::FieldMatrix<DT,n,n> permloc;
		Dune::FieldMatrix<DT,n,n> permlocWell;
		RT p0;
  };

}
#endif
