#ifndef TESTPROBLEM_UPSS_HH
#define TESTPROBLEM_UPSS_HH

#include "dumux/fractionalflow/fractionalflowmultiscaleproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class Grid, class Scalar, class VC, class CoarseScaleParameters> class UpsSProblem: public FractionalFlowMSProblem<
        Grid, Scalar, VC,CoarseScaleParameters>
{

    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    UpsSProblem(VC& variables, Fluid& wettingphase, Fluid& nonwettingphase, Matrix2p<Grid, Scalar>& soil, CoarseScaleParameters& coarseParams, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>), const FieldVector<Scalar,dim> Left = 0,
            const FieldVector<Scalar,dim> Right = 0, const bool capillarity = false) :
    FractionalFlowMSProblem<Grid, Scalar, VC, CoarseScaleParameters>(variables, wettingphase, nonwettingphase, soil, coarseParams, materialLaw, capillarity), Left_(Left[0]),
    Right_(Right[0]), eps_(1e-8)
    {}

    virtual Scalar sourcePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos)
    {
        return 0;
    }

    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos,
            const Element& element, const LocalPosition& localPos) const
    {
        if ((globalPos[0] < eps_ ) || globalPos[0]> (Right_ - eps_))
        return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_ )
        return Dune::BoundaryConditions::dirichlet;
        else
        return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_)
        return 1;
        // all other boundaries
        return 0;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_)
        return 1;
        // all other boundaries
        return 0;
    }

    Scalar neumannPress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        //        if (globalPos[0] > Right_ - eps_)
        //            return 3e-7;
        return 0;
    }

    Scalar neumannSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos, Scalar factor) const
    {
        if (globalPos[0]> (Right_ - eps_))
        return factor;
        // all other boundaries
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        //      if (globalPos[0] < eps_)
        //    return 1;
        //      else
        return 0;
    }

private:
    Scalar Left_;
    Scalar Right_;

    Scalar eps_;
};
}

#endif
