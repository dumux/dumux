// $Id$

#ifndef UNIFORMPROBLEM_HH
#define UNIFORMPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"
#include<dumux/material/twophaserelations.hh>
#include <dumux/material/property_baseclasses.hh>
#include <dune/istl/bvector.hh>

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class Grid, class Scalar, class VC>
class UniformProblem: public DiffusionProblem<Grid, Scalar, VC>
{
protected:
    enum
    {
        dim = Grid::dimension, dimWorld = Grid::dimensionworld
    };
typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    UniformProblem(VC& variables,Fluid& wettingPhase, Fluid& nonwettingPhase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>), const bool capillarity = false)
    : DiffusionProblem<Grid,Scalar,VC>(variables, wettingPhase, nonwettingPhase, soil, materialLaw, capillarity)
    {}

    UniformProblem()
    : DiffusionProblem<Grid,Scalar,VC>()
    {}

    Scalar sourcePress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos)
    {
        return 0;
    }

    typename Dune::BoundaryConditions::Flags bctypePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0]> 1-1E-6 || globalPos[0] < 1e-6)
        return Dune::BoundaryConditions::dirichlet;
        // all other boundaries
        return Dune::BoundaryConditions::neumann;
    }

    Scalar dirichletPress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return (globalPos[0] < 1e-6) ? 2e5 : 1e5;
    }

    Scalar dirichletSat(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] < 1e-6)
        return 0.8;
        else
        return 0.2;
    }

    Scalar neumannPress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        return 0;
    }
};
}

#endif
