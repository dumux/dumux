// $Id: lensproblem.hh 566 2008-09-11 11:38:31Z bernd $

#ifndef DUNE_LENSPROBLEM_HH
#define DUNE_LENSPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/property_baseclasses.hh>
#include<dumux/twophase/twophaseproblem.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the TwoPhase problem
 * @author Bernd Flemisch, Klaus Mosthaf
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
  class LensProblem : public TwoPhaseProblem<G, RT> {
	typedef typename G::ctype DT;
	enum {dim=G::dimension, m=2};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

  public:
	// Constructor
	LensProblem(Fluid& liq1, Fluid& liq2, Matrix2p<G, RT>& soil,
			const FieldVector<DT,dim>& outerLowerLeft = 0., const FieldVector<DT,dim>& outerUpperRight = 0,
			const FieldVector<DT,dim>& innerLowerLeft = 0., const FieldVector<DT,dim>& innerUpperRight = 0,
			TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G, RT>))
	: TwoPhaseProblem<G,RT>(liq1, liq2, soil, law),
		outerLowerLeft_(outerLowerLeft), outerUpperRight_(outerUpperRight),
		innerLowerLeft_(innerLowerLeft), innerUpperRight_(innerUpperRight),
		eps_(1e-8*outerUpperRight[0]),
		densityW_(liq1.density()), densityN_(liq2.density())
	{
		height_ = outerUpperRight[1] - outerLowerLeft[1];
		width_ = outerUpperRight[0] - outerLowerLeft[0];

		gravity_[0] = 0;
		gravity_[1] = -9.81;
	}

    // indices for the two phases
	enum {pWIdx = 0, sNIdx = 1};

	// function returning the BOUNDARY CONDITION TYPE depending on the position
	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,dim>& x, const Entity& e,
					const IntersectionIterator& intersectionIt,
					   const FieldVector<DT,dim>& xi) const
	{
		// boundary condition type is set to neumann
		FieldVector<BoundaryConditions::Flags, m> values(BoundaryConditions::neumann);

		if (x[0] < outerLowerLeft_[0] + eps_ || x[0] > outerUpperRight_[0] - eps_)
		{
			values = BoundaryConditions::dirichlet;
		}

		return values;
	}

	// definition of DIRICHLET boundary conditions depending on the position
	virtual FieldVector<RT,m> g (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		if (x[0] < outerLowerLeft_[0] + eps_)
		{
			RT a = -(1 + 0.5/height_);
			RT b = -a*outerUpperRight_[1];
			values[pWIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
			values[sNIdx] = 0;
		}
		else {
			RT a = -1;
			RT b = outerUpperRight_[1];
			values[pWIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
			values[sNIdx] = 0;
		}

		return values;
	}

	// definition of NEUMANN boundary conditions depending on the position
	virtual FieldVector<RT,m> J (const FieldVector<DT,dim>& x, const Entity& e,
				const IntersectionIterator& intersectionIt,
				  const FieldVector<DT,dim>& xi) const
	{
		FieldVector<RT,m> values(0);

		RT lambda = (outerUpperRight_[0] - x[0])/width_;
		if (lambda > 0.5 && lambda < 2.0/3.0 && x[1] > outerUpperRight_[1] - eps_) {
			values[sNIdx] = -0.04;
		}

		return values;
	}

	// definition of INITIAL VALUES for pressure and saturation
	virtual FieldVector<RT,m> initial (const FieldVector<DT,dim>& x, const Entity& e,
				  const FieldVector<DT,dim>& xi) const
	{

		FieldVector<RT,m> values;

		values[pWIdx] = -densityW_*gravity_[1]*(height_ - x[1]);

		if (x[0] < outerLowerLeft_[0] + eps_) {
			RT a = -(1 + 0.5/height_);
			RT b = -a*outerUpperRight_[1];
			values[pWIdx] = -densityW_*gravity_[1]*(a*x[1] + b);
		}

		if (x[0] > innerLowerLeft_[0] && x[0] < innerUpperRight_[0]
		    && x[1] > innerLowerLeft_[1] && x[1] < innerUpperRight_[1])
			values[sNIdx] = 0.0;//innerSnr_;
		else
			values[sNIdx] = 0.0;//outerSnr_;

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
		FieldVector<DT,dim> outerLowerLeft_;
		FieldVector<DT,dim> outerUpperRight_;
		FieldVector<DT,dim> innerLowerLeft_;
		FieldVector<DT,dim> innerUpperRight_;
		DT width_, height_;
		DT eps_;
		RT densityW_, densityN_;
		FieldVector<DT,dim> gravity_;
  };

} // end namespace
#endif
