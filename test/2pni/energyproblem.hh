// $Id: energyproblem.hh 882 2008-12-04 09:05:55Z melanie $

#ifndef DUNE_TWOPHEATPROBLEM_HH
#define DUNE_TWOPHEATPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/property_baseclasses.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/relperm_pc_law.hh>
#include<dumux/2pni/2pniproblem.hh>


/**
 * @file
 * @brief  Class for defining a nonisothermal two-phase problem
 * @author Bernd Flemisch, Klaus Mosthaf, Melanie Darcis
 */

namespace Dune
{
  template<class G, class RT>
  class TwoPHeatProblem : public TwoPhaseHeatProblem<G, RT> {
	typedef typename G::ctype DT;
	enum {dim=G::dimension, m=3};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
	// Constructor
	TwoPHeatProblem(Fluid& liq1, Fluid& liq2, Matrix2p<G, RT>& soil,
			TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G, RT>))
	: TwoPhaseHeatProblem<G,RT>(liq1, liq2, soil, law),
		eps_(1e-4)
	{
		gravity_[0] = 0;
		gravity_[1] = -9.81;
	}

    // indices for the two phases
	enum {pWIdx = 0, sNIdx = 1, teIdx = 2};

	// function returning the BOUNDARY CONDITION TYPE depending on the position
	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,dim>& x, const Entity& e,
					const IntersectionIterator& intersectionIt,
					   const FieldVector<DT,dim>& xi) const
	{
		// boundary condition type is set to neumann
		FieldVector<BoundaryConditions::Flags, m> values(BoundaryConditions::neumann);

		if (x[0] < eps_)
		{
			values = BoundaryConditions::dirichlet;
		}
		if (x[1] < 5 && x[0] > 10 - eps_)
		{
			values[teIdx] = BoundaryConditions::dirichlet;
		}
		return values;
	}

	// definition of DIRICHLET boundary conditions depending on the position
	virtual FieldVector<RT,m> g (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

			values[pWIdx] = 101300 - gravity_[1] * 1066.7 * (800 - x[1]);
			values[sNIdx] = 0;
			values[teIdx] = 283.15 + (800 - x[1])*0.03;

		return values;
	}

	// definition of NEUMANN boundary conditions depending on the position
	virtual FieldVector<RT,m> J (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		if (x[1] < 5 && x[0] > 10 - eps_) {
			values[pWIdx] = 0.0;
			values[sNIdx] = -0.04;
			values[teIdx] = 0.0;
		}

		return values;
	}

	// definition of INITIAL VALUES for pressure and saturation
	virtual FieldVector<RT,m> initial (const FieldVector<DT,dim>& x, const Entity& e,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values;

		values[pWIdx] = 101300 - gravity_[1] * 1066.7 * (800 - x[1]);
		values[sNIdx] = 0;
		values[teIdx] = 283.15 + (800 - x[1])*0.03;

		return values;
	}


	// function returning SOURCE/SINK terms
	virtual FieldVector<RT,m> q (const FieldVector<DT,dim>& x, const Entity& e,
					const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		return values;
	}

	virtual FieldVector<RT,dim> gravity () const
	{
		return gravity_;
	}

	private:
		DT eps_;
		FieldVector<DT,dim> gravity_;
  };

} // end namespace
#endif
