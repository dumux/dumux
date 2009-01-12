// $Id$
#ifndef DUNE_CO2PROBLEM_HH
#define DUNE_CO2PROBLEM_HH

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
   *    Template parameters are:
   *
   *    - Grid  a DUNE grid type
   *    - RT    type used for return values
   */
  template<class G, class RT>
  class CO2Problem : public TwoPhaseProblem<G, RT> {
    typedef typename G::ctype DT;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;
    enum {dim=G::dimension, m=2};

  public:
    enum {pWIdx = 0, sNIdx = 1};

    virtual FieldVector<RT,m> q (const FieldVector<DT,dim>& x, const Entity& e,
                    const FieldVector<DT,dim>& xi) const
    {
        FieldVector<RT,m> values(0);

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,dim>& x, const Entity& e,
                    const IntersectionIterator& intersectionIt,
                       const FieldVector<DT,dim>& xi) const
    {
        FieldVector<BoundaryConditions::Flags, m> values(Dune::BoundaryConditions::neumann);

        if(x[0] >= 300-1e-2 )
        values = Dune::BoundaryConditions::dirichlet;

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

        values[0] = 1.013e5 + (depthBOR_ - x[1]) * 1045 * 9.81;
        values[1] = 0.0;

        return values;
    }

    virtual FieldVector<RT,m> J (const FieldVector<DT,dim>& x, const Entity& e,
                const IntersectionIterator& intersectionIt, const FieldVector<DT,dim>& xi) const
    {
        FieldVector<RT,m> values(0);
        if(x[0] < 1.e-2 && x[1] < 30.)
        {
            values[0] = 0.0;//-4.046e-5;
            values[1] = -0.02;
        }
        return values;
    }

    // Initial Conditions for global vector x, element e and local vector xi
    virtual FieldVector<RT,m> initial (const FieldVector<DT,dim>& x, const Entity& e,
                  const FieldVector<DT,dim>& xi) const
    {
        FieldVector<RT,m> values(0);
        values[0] = 1.013e5 + (depthBOR_ - x[1]) * 1045 * 9.81;
        values[1] = 0.00;


        return values;
    }

    FieldVector<RT,dim> gravity () const
    {
        FieldVector<RT,dim> values(0);

        values[1] = -9.81;

        return values;
    }

    CO2Problem(Fluid& liq1, Fluid& liq2, Matrix2p<G, RT>& soil,
            TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G, RT>), RT depthBOR = 0.0)
    : TwoPhaseProblem<G,RT>(liq1, liq2, soil, law)

    {
            depthBOR_ = depthBOR;
    }

    private:
        RT depthBOR_;
  };

} // end namespace
#endif
