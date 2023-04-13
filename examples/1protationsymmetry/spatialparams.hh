// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_ONEP_ROTATION_SYMMETRY_SPATIAL_PARAMS_HH
#define DUMUX_ONEP_ROTATION_SYMMETRY_SPATIAL_PARAMS_HH

// ## Parameter distributions (`spatialparams.hh`)
//
// This file contains the __spatial parameters class__ which defines the
// distributions for the porous medium parameters permeability and porosity
// over the computational grid.
// [[content]]
// We include the spatial parameters class for single-phase models discretized
// by finite volume schemes, from which the spatial parameters defined for this
// example inherit.
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

// ### The spatial parameters class
//
// In the `RotSymExampleSpatialParams` class, we define the functions needed to describe
// the porous medium, that is, porosity and permeability.
// We inherit from the `FVPorousMediumFlowSpatialParamsOneP` class here, which is the base class
// for spatial parameters in the context of single-phase porous medium flow
// applications using finite volume discretization schemes.
// [[codeblock]]
namespace Dumux {

template<class GridGeometry, class Scalar>
class RotSymExampleSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, RotSymExampleSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = RotSymExampleSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
public:
    // Spatial parameter classes for porous medium flow applications need to
    // export the type used for intrinsic permeabilities.
    using PermeabilityType = Scalar;

    // In the constructor we obtain the permeability value from the input file.
    RotSymExampleSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    { permeability_ = getParam<Scalar>("SpatialParams.Permeability"); }
    // [[/codeblock]]

    // #### Porosity distribution
    // This function is used to define the porosity distribution in the
    // computational domain. Here, we use a constant porosity of 1.0.
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    // #### Permeability distribution
    // This function is used to define the permeability distribution in the
    // computational domain. Here, we use a constant permeability that is
    // defined in the input file.
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    { return permeability_; }

private:
    Scalar permeability_;
};

} // end namespace Dumux
// [[/content]]
#endif
