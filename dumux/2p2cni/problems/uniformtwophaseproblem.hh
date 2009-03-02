// $Id$

#ifndef DUNE_UNIFORMTWOPHASEPROBLEM_HH
#define DUNE_UNIFORMTWOPHASEPROBLEM_HH

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
#include<dumux/twophase/twophaseproblem_deprecated.hh>

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
class UniformTwoPhaseProblem : public TwoPhaseProblem<G, RT> {
    typedef typename G::ctype DT;
    enum {n=G::dimension, m=2};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;

public:
    virtual const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e,
                                          const FieldVector<DT,n>& xi)
    {
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

        if (x[0] > 600-1E-6)// || x[0] < 1e-6)
            //if (x[0] < 1E-6)
            values = Dune::BoundaryConditions::dirichlet;

        return values;
    }

    virtual FieldVector<RT,m> g (const FieldVector<DT,n>& x, const Entity& e,
                                 const IntersectionIterator& intersectionIt,
                                 const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,m> values(0);
        if (x[0] < 1e-6) {
            values[0] = 1e6;
            values[1] = 1;
        }
        else {
            values[0] = 0;
            values[1] = 0;
        }

        return values;
    }

    virtual FieldVector<RT,m> J (const FieldVector<DT,n>& x, const Entity& e,
                                 const IntersectionIterator& intersectionIt,
                                 const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,m> values(0);
        if (x[0] < 1e-6)
            values[0] = 1;

        return values;
    }

    virtual FieldVector<RT,m> initial (const FieldVector<DT,n>& x, const Entity& e,
                                       const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,m> values(0);
        values[0] = 1e6 - 1.0/600.0*1e6*x[0];
        //values[1] = 1 - 1.0/600.0*x[0];
        if (x[0] < 1e-6) {
            values[1] = 1e3;
        }

        return values;
    }

    double porosity (const FieldVector<DT,n>& x, const Entity& e,
                     const FieldVector<DT,n>& xi) const
    {
        return 1.0;
    }

    virtual FieldVector<RT,n> gravity () const
    {
        FieldVector<RT,n> values(0);

        return values;
    }

    virtual FieldVector<RT,4> materialLawParameters (const FieldVector<DT,n>& x, const Entity& e,
                                                     const FieldVector<DT,n>& xi) const
    {
        FieldVector<RT,4> values(0);

        return values;
    }

    UniformTwoPhaseProblem(TwoPhaseRelations& law = *(new DeprecatedLinearLaw))
        : TwoPhaseProblem<G, RT>(law)
    {
        permloc = 0;

        for (int i = 0; i < n; i++)
            permloc[i][i] = 1.0;
    }

private:
    Dune::FieldMatrix<DT,n,n> permloc;
};

}
#endif
