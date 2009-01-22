// $Id$
#ifndef CONVECTIONSUBPROBLEM_HH
#define CONVECTIONSUBPROBLEM_HH

#include "dumux/upscaledsaturation/preprocess/fractionalflowproblemsubprobs.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class Grid, class Scalar, class VC> class ConvSubProblemX1: public FractionalFlowProblemSubProbs<
        Grid, Scalar, VC>
{
    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

typedef    typename Grid::Traits::template Codim<0>::Entity Element;

public:
    ConvSubProblemX1(VC& variableobj, Fluid& wp, Fluid& nwp, Matrix2p<typename Grid::HostGridType, Scalar>& s, TwoPhaseRelations<typename Grid::HostGridType, Scalar>& law = *(new TwoPhaseRelations<typename Grid::HostGridType,Scalar>), GlobalPosition Left = 0,
            GlobalPosition Right = 0) :
    FractionalFlowProblemSubProbs<Grid, Scalar, VC>(variableobj, wp, nwp, s, law), LowerLeft_(Left),
    UpperRight_(Right), eps_(1e-6)
    {}

    virtual Scalar sourcePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos)
    {
        return 0;
    }

    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        if (globalPos[0] < LowerLeft_[0] + eps_ || globalPos[0]> (UpperRight_[0] - eps_))
        return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < LowerLeft_[0] + eps_)
        return Dune::BoundaryConditions::dirichlet;
        // all other boundaries
        return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < LowerLeft_[0] + eps_)
        return 1;
        // all other boundaries
        return 0;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < LowerLeft_[0] + eps_)
        return 1;
        // all other boundaries
        return 0;
    }

    Scalar neumannPress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0;
    }

    Scalar neumannSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos, Scalar& factor) const
    {
        if (globalPos[0]> (UpperRight_[0] - eps_))
        return factor;
        // all other boundaries
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0;
    }

private:
    GlobalPosition LowerLeft_;
    GlobalPosition UpperRight_;

    Scalar eps_;
};

template<class Grid, class Scalar, class VC> class ConvSubProblemY1: public FractionalFlowProblemSubProbs<
Grid , Scalar, VC>
{
    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename Grid::Traits::template Codim<0>::Entity Element;

public:
    ConvSubProblemY1(VC& variableobj, Fluid& wp, Fluid& nwp, Matrix2p<typename Grid::HostGridType, Scalar>& s, TwoPhaseRelations<typename Grid::HostGridType, Scalar>& law = *(new TwoPhaseRelations<typename Grid::HostGridType,Scalar>), GlobalPosition Left = 0,
            GlobalPosition Right = 0) :
    FractionalFlowProblemSubProbs<Grid, Scalar, VC>(variableobj, wp, nwp, s, law), LowerLeft_(Left),
    UpperRight_(Right), eps_(1e-8)
    {}

    virtual Scalar sourcePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos)
    {
        return 0;
    }

    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        if ((globalPos[1] < LowerLeft_[1] + eps_ || globalPos[1]> (UpperRight_[1] - eps_)))
        return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[1] < LowerLeft_[1] + eps_)
        return Dune::BoundaryConditions::dirichlet;
        else
        return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[1] < LowerLeft_[1] + eps_)
        return 1;
        // all other boundaries
        return 0;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[1] < LowerLeft_[1] + eps_)
        return 1;
        // all other boundaries
        return 0;
    }

    Scalar neumannPress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0;
    }

    Scalar neumannSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos, Scalar& factor) const
    {
        if (globalPos[1]> (UpperRight_[1] - eps_))
        return factor;
        // all other boundaries
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0;
    }

private:
    GlobalPosition LowerLeft_;
    GlobalPosition UpperRight_;

    Scalar eps_;
};
template<class Grid, class Scalar, class VC> class ConvSubProblemX2: public FractionalFlowProblemSubProbs<
Grid, Scalar, VC>
{
    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename Grid::Traits::template Codim<0>::Entity Element;

public:
    ConvSubProblemX2(VC& variableobj, Fluid& wp, Fluid& nwp, Matrix2p<typename Grid::HostGridType, Scalar>& s, TwoPhaseRelations<typename Grid::HostGridType, Scalar>& law = *(new TwoPhaseRelations<typename Grid::HostGridType,Scalar>), GlobalPosition Left = 0,
            GlobalPosition Right = 0) :
    FractionalFlowProblemSubProbs<Grid, Scalar, VC>(variableobj, wp, nwp, s, law), LowerLeft_(Left),
    UpperRight_(Right), eps_(1e-6)
    {}

    virtual Scalar sourcePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos)
    {
        return 0;
    }

    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        if (globalPos[0] < LowerLeft_[0] + eps_ || globalPos[0]> (UpperRight_[0] - eps_))
        return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0]> UpperRight_[0] - eps_)
        return Dune::BoundaryConditions::dirichlet;
        // all other boundaries
        return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0]> UpperRight_[0] - eps_)
        return 1;
        // all other boundaries
        return 0;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0]> UpperRight_[0] - eps_)
        return 1;
        // all other boundaries
        return 0;
    }

    Scalar neumannPress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0;
    }

    Scalar neumannSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos, Scalar& factor) const
    {
        if (globalPos[0] < (LowerLeft_[0] + eps_))
        return factor;
        // all other boundaries
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0;
    }

private:
    GlobalPosition LowerLeft_;
    GlobalPosition UpperRight_;

    Scalar eps_;
};

template<class Grid, class Scalar, class VC> class ConvSubProblemY2: public FractionalFlowProblemSubProbs<
Grid , Scalar, VC>
{
    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    typedef typename Grid::Traits::template Codim<0>::Entity Element;

public:
    ConvSubProblemY2(VC& variableobj, Fluid& wp, Fluid& nwp, Matrix2p<typename Grid::HostGridType, Scalar>& s, TwoPhaseRelations<typename Grid::HostGridType, Scalar>& law = *(new TwoPhaseRelations<typename Grid::HostGridType,Scalar>), GlobalPosition Left = 0,
            GlobalPosition Right = 0) :
    FractionalFlowProblemSubProbs<Grid, Scalar, VC>(variableobj, wp, nwp, s, law), LowerLeft_(Left),
    UpperRight_(Right), eps_(1e-8)
    {}

    virtual Scalar sourcePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos)
    {
        return 0;
    }

    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        if ((globalPos[1] < LowerLeft_[1] + eps_ || globalPos[1]> (UpperRight_[1] - eps_)))
        return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[1]> UpperRight_[1] - eps_)
        return Dune::BoundaryConditions::dirichlet;
        else
        return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[1]> UpperRight_[1] - eps_)
        return 1;
        // all other boundaries
        return 0;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[1]> UpperRight_[1] - eps_)
        return 1;
        // all other boundaries
        return 0;
    }

    Scalar neumannPress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0;
    }

    Scalar neumannSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos, Scalar& factor) const
    {
        if (globalPos[1] < (LowerLeft_[1] + eps_))
        return factor;
        // all other boundaries
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0;
    }

private:
    GlobalPosition LowerLeft_;
    GlobalPosition UpperRight_;

    Scalar eps_;
};
}
#endif
