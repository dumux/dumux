// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_EXAMPLE_SHALLOWWATER_FRICTION_PROPERTIES_HH
#define DUMUX_EXAMPLE_SHALLOWWATER_FRICTION_PROPERTIES_HH

// ## Compile-time settings (`properties.hh`)
//
// In this file, the type tag used for the shallow water flow simulation is defined,
// for which we then specialize `properties` to the needs of the desired setup.
//
// [[content]]
//
// ### Includes
// [[details]] include files
//
// The `ShallowWater` type tag specializes most of the `properties` required for a
// shallow water flow simulation in DuMu<sup>x</sup>. We will use this in the following to inherit the
// respective properties and subsequently specialize those `properties` for our
// type tag, which we want to modify or for which no meaningful default can be set.
#include <dumux/freeflow/shallowwater/model.hh>

// We want to use `YaspGrid`, an implementation of the dune grid interface for structured grids:
#include <dune/grid/yaspgrid.hh>

// In this example, we want to discretize the equations with the cell centered finite volume
// scheme using two-point-flux approximation:
#include <dumux/discretization/cctpfa.hh>

// We include the problem and spatial parameters headers used for this simulation.
#include "problem.hh"
#include "spatialparams.hh"
// [[/details]]
//
// ### Type tag definition
//
// First, a so-called type tag is created. Properties are traits specialized for this type tag (a simple `struct`).
// The properties of two other type tags are inherited by adding the alias `InheritsFrom`.
// Here, properties from the shallow water model (`TTag::ShallowWater`) and the
// cell-centered finite volume scheme with two-point-flux approximation (`TTag::CCTpfaModel`)
// are inherited. These other type tag definitions can be found in the included
// headers `dumux/freeflow/shallowwater/model.hh` and `dumux/discretization/cctpfa.hh`.
// [[codeblock]]
namespace Dumux::Properties {
namespace TTag {
struct RoughChannel { using InheritsFrom = std::tuple<ShallowWater, CCTpfaModel>; };
} // end namespace TTag
// [[/codeblock]]

// ### Property specializations
//
// We use a structured Cartesian grid with tensor product structure.
// `Dune::YaspGrid` (Yet Another Structure Parallel Grid) is defined in `dune/grid/yaspgrid.hh`
// in the Dune module `dune-grid`.
template<class TypeTag>
struct Grid<TypeTag, TTag::RoughChannel>
{ using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Next, we specialize the properties `Problem` and `SpatialParams` for our new type tag and
// set the type to our problem and spatial parameter classes implemented
// in `problem.hh` and `spatialparams.hh`.
// [[codeblock]]
template<class TypeTag>
struct Problem<TypeTag, TTag::RoughChannel>
{ using type = Dumux::RoughChannelProblem<TypeTag>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::RoughChannel>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;

    using type = RoughChannelSpatialParams<GridGeometry, Scalar, VolumeVariables>;
};
// [[/codeblock]]

// Finally, we enable caching for the grid geometry. When this feature
// is enabled, the entire finite-volume grid is precomputed and stored
// instead of preparing element-local geometries on the fly when assembling
// the linear system. This speeds up the simulation at the cost of a larger
// memory footprint.
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RoughChannel>
{ static constexpr bool value = true; };

} // end namespace Dumux::Properties
// [[/content]]
#endif
