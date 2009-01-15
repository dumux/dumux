// $Id$
#ifndef DUNE_TISSUE_TUMOR_PROBLEM_HH
#define DUNE_TISSUE_TUMOR_PROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>

#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/phaseproperties/phaseproperties1p.hh>
#include<dumux/1p2c/1p2cproblem.hh>
#include"tissue_soilproperties.hh"

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
  class TissueTumorProblem : public OnePTwoCProblem<G, RT> {
	typedef typename G::ctype DT;
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;
	enum {dim=G::dimension, m=2};
	enum {konti = 0, transport = 1};	// Solution vector index

  public:

	virtual FieldVector<RT,m> q (const FieldVector<DT,dim>& x, const Entity& e,
					const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);
		
		if(x[0]>10 && x[0]<12 && x[1]>10 && x[1]<12)
			values[0]=1.5e-6;
		return values;
	}

	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,dim>& x, const Entity& e,
					const IntersectionIterator& intersectionIt,
					   const FieldVector<DT,dim>& xi) const
	{
		FieldVector<BoundaryConditions::Flags, m> values(Dune::BoundaryConditions::dirichlet);

//		if( x[0]<1e-1 || x[0]> 22-1e-1 )
//			values = Dune::BoundaryConditions::dirichlet;
		
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

	virtual FieldVector<RT,m> g (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		if(x[0]<1e-1)
		 {
		  values[0] = -931;
		  values[1] = 1.1249e-8;
		 }
		  		
		else if(x[0]>22-1e-1)
         {
		 values[0] = -1067;  	//Dirichlet RB für Druck
		 values[1] = 1.1249e-8;       //Dirichlet RB für mol fraction x
		 }
			
		else
		{
		values[0] =-931-((136/22)*x[0]); //AB für Druck p
		values[1] = 1.1249e-8; //AB für mole fraction x	
		}
		return values;
	}

	virtual FieldVector<RT,m> J (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt, const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);
//		if(x[0] < 1.e-1)
//		{
//			values[0] = -3.8676e-8;
//			values[1] = -4.35064e-16;
//		}
		return values;
	}

	// Initial Conditions for global vector x, element e and local vector xi
	virtual FieldVector<RT,m> initial (const FieldVector<DT,dim>& x, const Entity& e,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);
				
		values[0] =-931-((136/22)*x[0]); //AB für Druck p
		
//		values[0] = -1067;
//		values[1] = 1.1249e-8;
		if(x[0]>10 && x[0]<12 && x[1]>10 && x[1]<12)
			values[1]=0;	
		else
			values[1]=1.1249e-8;
		
		return values;
	}



	TissueTumorProblem(Liquid_GL& phase, Matrix2p<G, RT>& soil)
	: OnePTwoCProblem<G, RT>(phase, soil)
	{
		
	}	
  };

}
#endif
