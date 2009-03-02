#ifndef BENCHMARK3PROBLEM_HH
#define BENCHMARK3PROBLEM_HH

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
#include<dumux/material/multicomponentrelations.hh>
#include<dumux/material/relperm_pc_law.hh>
#include<dumux/2p2cni/2p2cniproblem.hh>


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
template<class Grid, class RT>
class Benchmark3Problem : public TwoPTwoCNIProblem<Grid, RT> {
    typedef typename Grid::ctype Scalar;
    enum {dim=Grid::dimension, m=3};
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

public:

    virtual FieldVector<RT,m> q (const FieldVector<Scalar,dim>& x, const Entity& e,
                                 const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<RT,m> values(0);

        RT m;
        m =2.573325096E-5;
        if (x[0] > 5399. && x[0] < 5508. && x[1] > 3285. && x[1] < 3393. && x[2] < 286.)
            values[1] = m;

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<Scalar,dim>& x, const Entity& e,
                                                              const IntersectionIterator& intersectionIt,
                                                              const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<BoundaryConditions::Flags, m> values(Dune::BoundaryConditions::dirichlet);

        if(intersectionIt->boundaryId() == 1){
            values = Dune::BoundaryConditions::neumann;
        }

        return values;

    }

    virtual void dirichletIndex(const FieldVector<Scalar,dim>& x, const Entity& e,
                                const IntersectionIterator& intersectionIt,
                                const FieldVector<Scalar,dim>& xi, FieldVector<int,m>& dirichletIndex) const
    {
        for (int i = 0; i < m; i++)
            dirichletIndex[i]=i;
        return;
    }

    virtual FieldVector<RT,m> g (const FieldVector<Scalar,dim>& x, const Entity& e,
                                 const IntersectionIterator& intersectionIt,
                                 const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<RT,m> values(0);
        values[0] = 1.013e5 + (depthBOR_ - x[2]) * 1037.24 * 9.81;
        values[1] = 0.0;
        values[2] = 283.15 + (depthBOR_ - x[2])*0.03;

        return values;
    }

    virtual FieldVector<RT,m> J (const FieldVector<Scalar,dim>& x, const Entity& e,
                                 const IntersectionIterator& intersectionIt, const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<RT,m> values(0);

        return values;
    }

    // Initial Conditions for global vector x, element e and local vector xi
    virtual FieldVector<RT,m> initial (const FieldVector<Scalar,dim>& x, const Entity& e,
                                       const FieldVector<Scalar,dim>& xi) const
    {
        FieldVector<RT,m> values(0);
        values[0] = 1.013e5 + (depthBOR_ - x[2]) * 1037.24 * 9.81;
        values[1] = 0.0;
        values[2] = 283.15 + (depthBOR_ - x[2])*0.03;

        return values;
    }

    int initialPhaseState (const FieldVector<Scalar,dim>& x, const Entity& e,
                           const FieldVector<Scalar,dim>& xi) const
    {

        enum {gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase states
        int state;

        state = waterPhase;

        return state;
    }

    virtual FieldVector<RT,dim> gravity () const
    {
        FieldVector<RT,dim> values(0);

        values[2] = -9.81;

        return values;
    }

    Benchmark3Problem(Liquid_GL& liq, Gas_GL& gas, Matrix2p<Grid, RT>& soil,
                      TwoPhaseRelations<Grid, RT>& law = *(new TwoPhaseRelations<Grid, RT>),
                      MultiComp& multicomp = *(new CWaterAir), RT depthBOR = 0.0)
        : TwoPTwoCNIProblem<Grid, RT>(liq, gas, soil, multicomp, law), depthBOR_(depthBOR)
    {
    }

private:
    RT depthBOR_;
};

}
#endif
