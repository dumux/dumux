// $Id$
#ifndef INJECTION_PROBLEM_HH
#define INJECTION_PROBLEM_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"

namespace Dune {
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class Grid, class Scalar, class VC> class InjectionProblem :
        public FractionalFlowProblem<Grid,Scalar,VC> {

    enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;

public:
    InjectionProblem(VC& variables, Fluid& wettingphase, Fluid& nonwettingphase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>), const GlobalPosition LowerLeft = 0,
                const GlobalPosition UpperRight = 0) :
        FractionalFlowProblem<Grid, Scalar, VC>(variables, wettingphase, nonwettingphase, soil, materialLaw), Left_(LowerLeft[0]),
        Right_(UpperRight[0]), eps_(1e-8), depthBOR_(1000)
    {
        gravity_[0] = 0;
        gravity_[1] = -9.81;
    }

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
        if (globalPos[0] < eps_)
            return Dune::BoundaryConditions::dirichlet;
        else
            return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletPress(const GlobalPosition& globalPos, const Element& element,
                          const LocalPosition& localPos) const {
        if (globalPos[0] < eps_)
            return (1.5 - this->wettingPhase.density() * gravity_[1] * (depthBOR_ - globalPos[1]));
        // all other boundaries
        return 0;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Element& element,
                        const LocalPosition& localPos) const {
        if (globalPos[0] < eps_)
            return 1;
        // all other boundaries
        return 0;
    }

    Scalar neumannPress(const GlobalPosition& globalPos, const Element& element,
                        const LocalPosition& localPos) const {
        if (globalPos[1] < 15 && globalPos[1] > 5)
            return -1.5e-3;
        return 0;
    }

    Scalar initSat (const GlobalPosition& globalPos, const Element& element,
                    const LocalPosition& localPos) const
    {
            return 1;
    }
    const FieldVector<Scalar,dimWorld>& gravity() const
    {
        return gravity_;
    }


private:
    FieldVector<Scalar,dimWorld> gravity_;

    Scalar Left_;
    Scalar Right_;
    Scalar eps_;
    Scalar depthBOR_;
};
}

#endif
