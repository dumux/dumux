#ifndef FIVESPOTDIFFPROBLEM_HH
#define FIVESPOTDIFFPROBLEM_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class Grid, class Scalar, class VC>
class FivespotProblemCase1: public FractionalFlowProblem<Grid, Scalar, VC>
{
    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
    FivespotProblemCase1(VC& variables, Fluid& wettingphase, Fluid& nonwettingphase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>),
            Scalar bcf = 11, const bool capillarity = false)
    : FractionalFlowProblem<Grid, Scalar, VC>(variables, wettingphase, nonwettingphase, soil, materialLaw, capillarity),
    LowerLeft_(variables.grid.lowerLeft()), UpperRight_(variables.grid.upperRight()),
    eps_(1e-8*UpperRight_[0]),bcf_(bcf),
    pressBC_(2e5)
    {}

    Scalar sourcePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos)
    {
        return 0;
    }

    typename BoundaryConditions::Flags bctypePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if ((globalPos[0] < LowerLeft_[0] + eps_ && globalPos[1] < LowerLeft_[1] + bcf_) || //lower left
                (globalPos[1] < LowerLeft_[1] + eps_ && globalPos[0] < LowerLeft_[0] + bcf_)) //lower left
        return BoundaryConditions::dirichlet;

        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if ((globalPos[0] < LowerLeft_[0] + eps_ && globalPos[1] < LowerLeft_[1] + bcf_) ||
                (globalPos[1] < LowerLeft_[1] + eps_ && globalPos[0] < LowerLeft_[0] + bcf_) ||
                (globalPos[0]> UpperRight_[0] - eps_ && globalPos[1]> UpperRight_[1] - bcf_) ||
                (globalPos[1]> UpperRight_[1] - eps_ && globalPos[0]> UpperRight_[0] - bcf_))
        return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    Scalar dirichletPress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if ((globalPos[0] < LowerLeft_[0] + eps_ && globalPos[1] < LowerLeft_[1] + bcf_) || //lower left
                (globalPos[1] < LowerLeft_[1] + eps_ && globalPos[0] < LowerLeft_[0] + bcf_)) //lower left
        return pressBC_;
        // all other boundaries
        return 0;
    }

    Scalar dirichletSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if ((globalPos[0] < LowerLeft_[0] + eps_ && globalPos[1] < LowerLeft_[1] + bcf_) ||
                (globalPos[1] < LowerLeft_[1] + eps_ && globalPos[0] < LowerLeft_[0] + bcf_))
        return 0.8;
        else
        return 0.2;
    }

    Scalar neumannPress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if ((globalPos[0]> UpperRight_[0] - eps_ && globalPos[1]> UpperRight_[1] - bcf_) ||
                (globalPos[1]> UpperRight_[1] - eps_ && globalPos[0]> UpperRight_[0] - bcf_))
        {
            if (bcf_ == 21)
            return 1e-6;//15x15 cells: bcf_ 21
            if (bcf_ == 11)
            return 2e-6; //30x30 cells: bcf_ 11
            if (bcf_ == 6)
            return 4e-6;// 60x60 cells: bcf_ 6
        }
        //all other boundaries
        return 0;

    }
    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0.2;
    }

private:
    FieldVector<Scalar,dim> LowerLeft_;
    FieldVector<Scalar,dim> UpperRight_;

    Scalar eps_;
    Scalar bcf_;
    Scalar pressBC_, prightbc_;
};
template<class Grid,class Scalar, class VC>
class FivespotProblemCase2 : public FractionalFlowProblem<Grid,Scalar,VC>
{
    enum
    {   dim=Grid::dimension, dimWorld = Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
    FivespotProblemCase2(VC& variables, Fluid& wettingphase, Fluid& nonwettingphase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>),
            Scalar bcf = 11, const bool capillarity = false)
    : FractionalFlowProblem<Grid, Scalar, VC>(variables, wettingphase, nonwettingphase, soil, materialLaw, capillarity),
    LowerLeft_(variables.grid.lowerLeft()), UpperRight_(variables.grid.upperRight()),
    eps_(1e-8*UpperRight_[0]),bcf_(bcf),
    pressBC_(2e5)
    {}

    Scalar sourcePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos)
    {
        return 0;
    }

    typename BoundaryConditions::Flags bctypePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if ((globalPos[0] < LowerLeft_[0] + eps_ && globalPos[1] < LowerLeft_[1] + bcf_) || //lower left
                (globalPos[1] < LowerLeft_[1] + eps_ && globalPos[0] < LowerLeft_[0] + bcf_) || //lower left
                (globalPos[0]> UpperRight_[0] - eps_ && globalPos[1]> UpperRight_[1] - bcf_) || //upper right
                (globalPos[1]> UpperRight_[1] - eps_ && globalPos[0]> UpperRight_[0] - bcf_)) //upper right
        return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if ((globalPos[0] < LowerLeft_[0] + eps_ && globalPos[1] < LowerLeft_[1] + bcf_) || //lower left
                (globalPos[1] < LowerLeft_[1] + eps_ && globalPos[0] < LowerLeft_[0] + bcf_) || //lower left
                (globalPos[0]> UpperRight_[0] - eps_ && globalPos[1]> UpperRight_[1] - bcf_) || //upper right
                (globalPos[1]> UpperRight_[1] - eps_ && globalPos[0]> UpperRight_[0] - bcf_) || //upper right
                (globalPos[0] < LowerLeft_[0] + eps_ && globalPos[1]> UpperRight_[1] - bcf_) || //upper left
                (globalPos[1]> UpperRight_[1] - eps_ && globalPos[0] < LowerLeft_[0] + bcf_) || //upper left
                (globalPos[0]> UpperRight_[0] - eps_ && globalPos[1] < LowerLeft_[1] + bcf_) || //lower right
                (globalPos[1] < LowerLeft_[1] + eps_ && globalPos[0]> UpperRight_[0] - bcf_)) //lower right
        return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    Scalar dirichletPress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if ((globalPos[0] < LowerLeft_[0] + eps_ && globalPos[1] < LowerLeft_[1] + bcf_) || //lower left
                (globalPos[1] < LowerLeft_[1] + eps_ && globalPos[0] < LowerLeft_[0] + bcf_) || //lower left
                (globalPos[0]> UpperRight_[0] - eps_ && globalPos[1]> UpperRight_[1] - bcf_) || //upper right
                (globalPos[1]> UpperRight_[1] - eps_ && globalPos[0]> UpperRight_[0] - bcf_)) //upper right
        return pressBC_;
        // all other boundaries
        return 0;
    }

    Scalar dirichletSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if ((globalPos[0] < LowerLeft_[0] + eps_ && globalPos[1] < LowerLeft_[1] + bcf_) || //lower left
                (globalPos[1] < LowerLeft_[1] + eps_ && globalPos[0] < LowerLeft_[0] + bcf_) || //lower left
                (globalPos[0]> UpperRight_[0] - eps_ && globalPos[1]> UpperRight_[1] - bcf_) || //upper right
                (globalPos[1]> UpperRight_[1] - eps_ && globalPos[0]> UpperRight_[0] - bcf_)) //upper right
        return 0.8;
        // all other boundaries
        return 0.2;
    }

    Scalar neumannPress (const GlobalPosition& globalPos, const Element& element,
            const GlobalPosition& globalPosi) const
    {
        if ((globalPos[0] < LowerLeft_[0] + eps_ && globalPos[1]> UpperRight_[1] - bcf_) || //upper left
                (globalPos[1]> UpperRight_[1] - eps_ && globalPos[0] < LowerLeft_[0] + bcf_) || //upper left
                (globalPos[0]> UpperRight_[0] - eps_ && globalPos[1] < LowerLeft_[1] + bcf_) || //lower right
                (globalPos[1] < LowerLeft_[1] + eps_ && globalPos[0]> UpperRight_[0] - bcf_)) //lower right

        {
            if (bcf_ == 21)
            return 1e-6;//15x15 cells: bcf_ 21
            if (bcf_ == 11)
            return 2e-6; //30x30 cells: bcf_ 11
            if (bcf_ == 6)
            return 4e-6;// 60x60 cells: bcf_ 6
        }
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0.2;
    }

private:
    FieldVector<Scalar,dim> LowerLeft_;
    FieldVector<Scalar,dim> UpperRight_;

    Scalar eps_;
    Scalar bcf_;
    Scalar pressBC_, prightbc_;
};

}

#endif
