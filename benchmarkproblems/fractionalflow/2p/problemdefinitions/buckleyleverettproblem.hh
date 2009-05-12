#ifndef DUNE_BUCKLEYLEVERETTPROBLEM_HH
#define DUNE_BUCKLEYLEVERETTPROBLEM_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief buckley-leverett problem

template<class Grid, class Scalar, class VC>
class BuckleyLeverettProblem: public FractionalFlowProblem<Grid, Scalar, VC>
{
    template<int dim> struct ElementLayout
    {
        bool contains(Dune::GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

typedef    typename Grid::ctype DT;
    enum
    {   dim = Grid::dimension, dimWorld = Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Entity;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:

    virtual Scalar sourcePress (const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos)
    {
        return 0;
    }

    BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos,
            const Entity& element, const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_)
        return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos,
            const Entity& element, const LocalPosition& localPos) const
    {
        if (globalPos[0]> (right_ - eps_) || globalPos[0] < eps_)
        return Dune::BoundaryConditions::dirichlet;
        else
        return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_)
        return pleftbc_;
        // all other boundaries
        return prightbc_;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_)
        return 0.8;
        // all other boundaries
        return 0.2;
    }

    Scalar neumannPress(const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0]> right_ - eps_)
        return 3e-7;
        //if (globalPos[0] < eps_) return -3e-7;
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Entity& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_)
        return 0.8;
        else
        return 0.2;
    }

public:
    BuckleyLeverettProblem(VC& variables, Fluid& wettingphase, Fluid& nonwettingphase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>&), const GlobalPosition Left = 0,
            const GlobalPosition Right = 0, Scalar pleftbc=2e5, Scalar prightbc=1.999986e5, const bool capillarity = false) :
    FractionalFlowProblem<Grid, Scalar, VC>(variables, wettingphase, nonwettingphase, soil, materialLaw, capillarity),
    left_(Left[0]), right_(Right[0]), pleftbc_(pleftbc), prightbc_(prightbc), eps_(1e-8)
    {}

private:
    Scalar left_, right_;
    Scalar pleftbc_, prightbc_;
    Scalar eps_;

};
}
#endif
