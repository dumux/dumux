// $Id:$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/

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
template<class GridView, class Scalar, class VC>
class UniformProblem: public DiffusionProblem<GridView, Scalar, VC>
{
protected:
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    typedef typename GridView::Grid Grid;
typedef    typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

public:
    UniformProblem(VC& variables,Fluid& wettingPhase, Fluid& nonwettingPhase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>), FieldVector<Scalar, dim> upperRight = *(new FieldVector<Scalar, dim>(1)))
    : DiffusionProblem<GridView,Scalar,VC>(variables, wettingPhase, nonwettingPhase, soil, materialLaw), upperRight_(upperRight)
    {}

    UniformProblem()
    : DiffusionProblem<GridView,Scalar,VC>()
    {}

    Scalar sourcePress(const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos)
    {
        return 0;
    }

    typename Dune::BoundaryConditions::Flags bctypePress (const GlobalPosition& globalPos, const Element& element,
            const LocalPosition& localPos) const
    {
        if (globalPos[0] > upperRight_[0] - 1e-6 || globalPos[0] < 1e-6)
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
private:
    FieldVector<Scalar, dim> upperRight_;
};
}

#endif
