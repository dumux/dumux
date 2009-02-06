// $Id$
#ifndef DUNE_BRINECO2PROBLEM_HH
#define DUNE_BRINECO2PROBLEM_HH

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
 *    - Scalar    type used for return values
 */
template<class Grid, class Scalar>
class BrineCO2Problem: public TwoPTwoCNIProblem<Grid, Scalar>
{
typedef    typename Grid::ctype CoordScalar;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;
    enum
    {   dim=Grid::dimension, numEq=3};
    enum
    {   wComp = 0, nComp = 1, temp = 2};

public:

    virtual FieldVector<Scalar,numEq> q (const FieldVector<CoordScalar,dim>& globalPos, const Element& element,
            const FieldVector<CoordScalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const FieldVector<CoordScalar,dim>& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const FieldVector<CoordScalar,dim>& localPos) const
    {
        FieldVector<BoundaryConditions::Flags, numEq> values(Dune::BoundaryConditions::neumann);

        if(globalPos[0] >= 10-1e-2 )
        values = Dune::BoundaryConditions::dirichlet;
        if(globalPos[0] < 1e-2 && globalPos[1] < 5)
        values[temp] = Dune::BoundaryConditions::dirichlet;

        return values;
    }

    virtual void dirichletIndex(const FieldVector<CoordScalar,dim>& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const FieldVector<CoordScalar,dim>& localPos, FieldVector<int,numEq>& dirichletIndex) const
    {
        for (int i = 0; i < numEq; i++)
        dirichletIndex[i]=i;
        return;
    }

    virtual FieldVector<Scalar,numEq> g (const FieldVector<CoordScalar,dim>& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const FieldVector<CoordScalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        values[wComp] = 1.013e5 + (depthBOR_ - globalPos[1]) * 1045 * 9.81;
        values[nComp] = 0.0;
        values[temp] = 283.15 + (depthBOR_ - globalPos[1])*0.03;

        return values;
    }

    virtual FieldVector<Scalar,numEq> J (const FieldVector<CoordScalar,dim>& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt, const FieldVector<CoordScalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);
        if(globalPos[0] <= 1.e-2 && globalPos[1] <= 5.)
        {
            values[wComp] = 0.0;//-4.046e-5;
            values[nComp] = -0.0002;
        }
        return values;
    }

    // Initial Conditions for global vector globalPos, element element and local vector localPos
    virtual FieldVector<Scalar,numEq> initial (const FieldVector<CoordScalar,dim>& globalPos, const Element& element,
            const FieldVector<CoordScalar,dim>& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);
        values[wComp] = 1.013e5 + (depthBOR_ - globalPos[1]) * 1045 * 9.81;
        values[nComp] = 0.00;
        values[temp] = 283.15 + (depthBOR_ - globalPos[1])*0.03;

        return values;
    }

    int initialPhaseState (const FieldVector<CoordScalar,dim>& globalPos, const Element& element,
            const FieldVector<CoordScalar,dim>& localPos) const
    {
        enum
        {   gasPhase = 0, waterPhase = 1, bothPhases = 2}; // Phase states
        int state;

        state = waterPhase;

        return state;
    }

    FieldVector<Scalar,dim> gravity () const
    {
        FieldVector<Scalar,dim> values(0);

        values[dim] = -9.81;

        return values;
    }

    BrineCO2Problem(Liquid_GL& liq, Gas_GL& gas, Matrix2p<Grid, Scalar>& soil,
            TwoPhaseRelations<Grid, Scalar>& law = *(new TwoPhaseRelations<Grid, Scalar>),
            MultiComp& multicomp = *(new CWaterAir), Scalar depthBOR = 0.0)
    : TwoPTwoCNIProblem<Grid, Scalar>(liq, gas, soil, multicomp, law)
    {
        depthBOR_ = depthBOR;
    }

private:
    Scalar depthBOR_, soilDens_, soilHeatCp_, soilLDry_, soilLSw_;
};

}
#endif
