#ifndef DUNE_MCWHORTERPROBLEM_HH
#define DUNE_MCWHORTERPROBLEM_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief McWhorter transport problem

template<class GridView, class Scalar, class VC>
class McWhorterProblem: public FractionalFlowProblem<GridView, Scalar, VC>
{
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
typedef    typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Entity;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:

    virtual Scalar sourcePress (const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos)
    {
        return 0;
    }

    BoundaryConditions::Flags bctypePress (const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_)
        return BoundaryConditions::dirichlet;

        // all other boundaries
        else
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_)// || globalPos[0]> right_ - eps_)
        return Dune::BoundaryConditions::dirichlet;
        else
        return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletPress (const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos) const
    {
        return 2e5;
    }

    Scalar dirichletSat (const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_)
        return 1.0;
        else
        return 0.0;
    }

    Scalar neumannPress (const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos) const
    {
        return 0.0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos) const
    {
        return 0.0;
    }

    McWhorterProblem(VC& variables, Fluid& wettingphase, Fluid& nonwettingphase,
            Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>&),
            const GlobalPosition Right = 0)
    : FractionalFlowProblem<GridView, Scalar, VC>(variables, wettingphase, nonwettingphase, soil, materialLaw),
    right_(Right[0]), eps_(1e-8)
    {}

private:
    Scalar right_;
    Scalar eps_;
 };
}
#endif
