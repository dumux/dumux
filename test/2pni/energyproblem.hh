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
#include<dumux/auxiliary/basicdomain.hh>

/**
 * @file
 * @brief  Class for defining a nonisothermal two-phase problem
 * @author Bernd Flemisch, Klaus Mosthaf, Melanie Darcis
 */

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class TwoPHeatProblem: public TwoPhaseHeatProblem<Grid, Scalar>
{
typedef    typename Grid::ctype CoordScalar;
    enum
    {   dim=Grid::dimension, numEq=3};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

    typedef BasicDomain<Grid, Scalar> ParentType;
    // the domain traits of the domain
    typedef typename ParentType::DomainTraits DomainTraits;
    typedef typename DomainTraits::LocalPosition LocalPosition;
    typedef typename DomainTraits::GlobalPosition GlobalPosition;

public:
    // Constructor
    TwoPHeatProblem(Fluid& liq1, Fluid& liq2, Matrix2p<Grid, Scalar>& soil,
            TwoPhaseRelations<Grid, Scalar>& law = *(new TwoPhaseRelations<Grid, Scalar>))
    : TwoPhaseHeatProblem<Grid,Scalar>(liq1, liq2, soil, law),
    eps_(1e-4)
    {
        gravity_[0] = 0;
        gravity_[1] = -9.81;
    }

    // indices for the two phases
    enum
    {   wPhase = 0, nPhase = 1, heat = 2};

    // function returning the BOUNDARY CONDITION TYPE depending on the position
    virtual FieldVector<BoundaryConditions::Flags, numEq> bctype (const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const LocalPosition& localPos) const
    {
        // boundary condition type is set to neumann
        FieldVector<BoundaryConditions::Flags, numEq> values(BoundaryConditions::neumann);

        if (globalPos[0] < eps_)
        {
            values = BoundaryConditions::dirichlet;
        }
        if (globalPos[1] < 5 && globalPos[0]> 10 - eps_)
        {
            values[heat] = BoundaryConditions::dirichlet;
        }
        return values;
    }

    // definition of DIRICHLET boundary conditions depending on the position
    virtual FieldVector<Scalar,numEq> g (const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const LocalPosition& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        values[wPhase] = 101300 - gravity_[1] * 1066.7 * (800 - globalPos[1]);
        values[nPhase] = 0;
        values[heat] = 283.15 + (800 - globalPos[1])*0.03;

        return values;
    }

    // definition of NEUMANN boundary conditions depending on the position
    virtual FieldVector<Scalar,numEq> J (const GlobalPosition& globalPos, const Element& element,
            const IntersectionIterator& intersectionIt,
            const LocalPosition& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        if (globalPos[1] < 5 && globalPos[0]> 10 - eps_)
        {
            values[wPhase] = 0.0;
            values[nPhase] = -0.04;
            values[heat] = 0.0;
        }

        return values;
    }

    // definition of INITIAL VALUES for pressure and saturation
    virtual FieldVector<Scalar,numEq> initial (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        FieldVector<Scalar,numEq> values;

        values[wPhase] = 101300 - gravity_[1] * 1066.7 * (800 - globalPos[1]);
        values[nPhase] = 0;
        values[heat] = 283.15 + (800 - globalPos[1])*0.03;

        return values;
    }

    // function returning SOURCE/SINK terms
    virtual FieldVector<Scalar,numEq> q (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        FieldVector<Scalar,numEq> values(0);

        return values;
    }

    virtual FieldVector<Scalar,dim> gravity () const
    {
        return gravity_;
    }

private:
    CoordScalar eps_;
    FieldVector<CoordScalar,dim> gravity_;
};

} // end namespace
#endif
