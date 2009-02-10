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
#include<dumux/auxiliary/basicdomain.hh>

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
class CO2Problem: public TwoPhaseProblem<Grid, Scalar>
{
typedef    typename Grid::ctype CoordScalar;
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;
    typedef BasicDomain<Grid, Scalar> ParentType;
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits DomainTraits;
    typedef typename DomainTraits::LocalPosition LocalPosition;
    typedef typename DomainTraits::GlobalPosition GlobalPosition;

    enum
    {   dim=Grid::dimension, numEq=2};

public:
    enum
    {   wPhase = 0, nPhase = 1};

    virtual FieldVector<Scalar,numEq> q (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        return values;
    }

    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& isIt,
            const LocalPosition& localPos) const
    {
        FieldVector<BoundaryConditions::Flags, numEq> values(Dune::BoundaryConditions::neumann);

        if(globalPos[0] >= 300-1e-2 )
        values = Dune::BoundaryConditions::dirichlet;

        return values;
    }

    virtual void dirichletIndex(const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& isIt,
            const LocalPosition& localPos, FieldVector<int,numEq>& dirichletIndex) const
    {
        for (int equationNumber = 0; equationNumber < numEq; equationNumber++)
        dirichletIndex[equationNumber]=equationNumber;
        return;
    }

    virtual FieldVector<Scalar,numEq> g (const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& isIt,
            const LocalPosition& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        values[wPhase] = 1.013e5 + (depthBOR_ - globalPos[1]) * 1045 * 9.81;
        values[nPhase] = 0.0;

        return values;
    }

    virtual FieldVector<Scalar,numEq> J (const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& isIt, const LocalPosition& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);
        if(globalPos[0] < 1.e-2 && globalPos[1] < 30.)
        {
            values[wPhase] = 0.0;//-4.046e-5;
            values[nPhase] = -0.02;
        }
        return values;
    }

    // Initial Conditions for global vector globalPos, element element and local vector localPos
    virtual FieldVector<Scalar,numEq> initial (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);
        values[wPhase] = 1.013e5 + (depthBOR_ - globalPos[1]) * 1045 * 9.81;
        values[nPhase] = 0.00;

        return values;
    }

    FieldVector<Scalar,dim> gravity () const
    {
        FieldVector<Scalar,dim> values(0);

        values[dim-1] = -9.81;

        return values;
    }

    CO2Problem(Fluid& liq1, Fluid& liq2, Matrix2p<Grid, Scalar>& soil,
            TwoPhaseRelations<Grid, Scalar>& law = *(new TwoPhaseRelations<Grid, Scalar>), Scalar depthBOR = 0.0)
    : TwoPhaseProblem<Grid,Scalar>(liq1, liq2, soil, law)

    {
        depthBOR_ = depthBOR;
    }

private:
    Scalar depthBOR_;
};

} // end namespace
#endif
