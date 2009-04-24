// $Id$
#ifndef TUTORIALPROBLEM_COUPLED_HH
#define TUTORIALPROBLEM_COUPLED_HH

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

namespace Dune
{

/** \todo Please doc me! */

template<class G, class RT> class TutorialProblemCoupled /*@\label{tutorial-coupled:tutorialproblem}@*/
    : public TwoPhaseProblem<G, RT>
{

    typedef typename G::ctype DT;
    enum {n=G::dimension, m=2};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename G::LeafGridView::IntersectionIterator
    IntersectionIterator;
public:
    TutorialProblemCoupled(Fluid& wp, Fluid& nwp, Matrix2p<G, RT>& s,
                           TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G, RT>),
                           const FieldVector<DT,n> Left = 0, const FieldVector<DT,n> Right = 0)
        : TwoPhaseProblem<G,RT>(wp, nwp, s, law),
          Left_(Left[0]), Right_(Right[0]), eps_(1e-8)
    {}

    // function returning source/sink terms for both mass balance equations
    // depending on the position within the domain
    virtual FieldVector<RT,m> q (const FieldVector<DT,n>& x, const Entity& e,  /*@\label{tutorial-coupled:q}@*/
                                 const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,m> values(0);

        return values;
    }

    // function returning the boundary condition type for the solution of
    // both mass balance equations depending on the position within the domain
    virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,n>& x, const Entity& e,
                                                              const IntersectionIterator& intersectionIt,
                                                              const FieldVector<DT,n>& xi) const /*@\label{tutorial-coupled:bctype}@*/
    {
        // boundary condition type is set to neumann
        FieldVector<BoundaryConditions::Flags, m> values(BoundaryConditions::neumann);

        // regions where dirichlet boundary conditions are applied
        if (x[0] < eps_)
        {
            values = BoundaryConditions::dirichlet;
        }

        return values;
    }

    // function returning the Dirichlet boundary condition for the solution of
    // both mass balance equations depending on the position within the domain
    virtual FieldVector<RT,m>dirichlet(const FieldVector<DT,n>& x, const Entity& e,  /*@\label{tutorial-coupled:g}@*/
                                 const IntersectionIterator& intersectionIt,
                                 const FieldVector<DT,n>& xi) const
    {
        // values are set to zero
        FieldVector<RT,m> values(0);

        values[0] = 2.e5; // pressure value
        values[1] = 0.0;  // saturation value

        return values;
    }

    // function returning the Neumann boundary condition for the solution
    // of the coupled system of balance equations depending on the position within the domain
    virtual FieldVector<RT,m> neumann(const FieldVector<DT,n>& x, const Entity& e,
                                 const IntersectionIterator& intersectionIt,
                                 const FieldVector<DT,n>& xi) const /*@\label{tutorial-coupled:J}@*/
    {
        // values are set to zero
        FieldVector<RT,m> values(0);

        // region where value for the water flux is not zero
        if (x[0]> Right_ - eps_)
        {

            values[0] = 0.;
            values[1] = 3e-4;
        }
        return values;
    }

    // function returning the initial values for pressure and saturation
    // depending on the position within the domain
    virtual FieldVector<RT,m> initial (const FieldVector<DT,n>& x, const Entity& e,  /*@\label{tutorial-coupled:initial}@*/
                                       const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,m> values(0);

        values[0] = 2.e5;
        values[1] = 1.0;
        return values;
    }

    virtual FieldVector<RT,n> gravity () const
    {
        FieldVector<DT,n> gravity(0);
        return gravity;
    }

private:
    DT Left_;
    DT Right_;

    RT eps_;
};

} // end namespace
#endif
