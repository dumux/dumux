// $Id$
#ifndef TEST_FRACTIONALFLOW_TESTPROBLEM_HH
#define TEST_FRACTIONALFLOW_TESTPROBLEM_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"

namespace Dune {
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class Grid, class Scalar, class VC> class FractionalFlowTestProblem :
        public FractionalFlowProblem<Grid,Scalar,VC> {

    enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
    FractionalFlowTestProblem(VC& variables, Fluid& wettingphase, Fluid& nonwettingphase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>), const GlobalPosition LowerLeft = 0,
                const GlobalPosition UpperRight = 0, const bool capillarity = false) :
        FractionalFlowProblem<Grid, Scalar, VC>(variables, wettingphase, nonwettingphase, soil, materialLaw, capillarity), Left_(LowerLeft[0]),
        Right_(UpperRight[0]), eps_(1e-8)
    {}

    virtual Scalar sourcePress  (const GlobalPosition& globalPos, const Element& element,
                                 const LocalPosition& localPos)
    {
        return 0;
    }

    typename BoundaryConditions::Flags bctypePress(const GlobalPosition& globalPos,
                                                   const Element& element, const LocalPosition& localPos) const {
        if ((globalPos[0] < eps_))
            return BoundaryConditions::dirichlet;
        // all other boundaries
        return BoundaryConditions::neumann;
    }

    BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
                                         const LocalPosition& localPos) const
    {
        if (globalPos[0] > (Right_ - eps_) || globalPos[0] < eps_)
            return Dune::BoundaryConditions::dirichlet;
        else
            return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Element& element,
                          const LocalPosition& localPos) const {
        if (globalPos[0] < eps_)
            return 2e5;
        // all other boundaries
        return 1.999985e5;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Element& element,
                        const LocalPosition& localPos) const {
        if (globalPos[0] < eps_)
            return 0.8;
        // all other boundaries
        return 0;
    }

    Scalar neumannPress(const GlobalPosition& globalPos, const Element& element,
                        const LocalPosition& localPos) const {
        if (globalPos[0] > Right_ - eps_)
            return 3e-7;
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
                    const LocalPosition& localPos) const
    {
        if (globalPos[0] < eps_)
            return 0.8;
        else
            return 0.2;
    }

private:
    Scalar Left_;
    Scalar Right_;

    Scalar eps_;
};
}

#endif
