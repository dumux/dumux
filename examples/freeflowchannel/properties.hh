// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_EXAMPLES_FREEFLOW_CHANNEL_PROPERTIES_HH
#define DUMUX_EXAMPLES_FREEFLOW_CHANNEL_PROPERTIES_HH

// ## Compile-time settings (`properties.hh`)
//
// In this file, the type tag used for this simulation is defined,
// for which we then specialize properties (compile time options) to the needs of the desired setup.
//
// [[content]]
//
// ### Includes
// [[details]] includes
//
// The `NavierStokes` type tag specializes most of the properties required for Navier-
// Stokes single-phase flow simulations in DuMuX. We will use this in the following to inherit the
// respective properties and subsequently specialize those properties for our
// type tag, which we want to modify or for which no meaningful default can be set.
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/navierstokes/staggered/problem.hh>

// We want to use `YaspGrid`, an implementation of the dune grid interface for structured grids:
#include <dune/grid/yaspgrid.hh>
// In this example, we want to discretize the equations with the staggered-grid
// scheme which is so far the only available option for free-flow models in DuMux:
#include <dumux/discretization/staggered/freeflow/properties.hh>
// The fluid properties are specified in the following headers (we use a liquid with constant properties as the fluid phase):
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

// We include the problem header used for this simulation.
#include "problem.hh"
// [[/details]]
//
// ### Type tag definition
//
// We define a type tag for our simulation with the name `ChannelExample`
// and inherit the properties specialized for the type tags `NavierStokes` and `StaggeredFreeFlowModel`.
// This way, most of the properties required for Navier-Stokes single-phase flow simulations
// using the staggered-grid scheme are conveniently specialized for our new type tag.
// However, some properties depend on user choices and no meaningful default value can be set.
// Those properties will be addressed later in this file.
// Please note that, in this example, we actually want to solve the Stokes instead of the
// Navier-Stokes equations. This can be achieved at runtime by setting the parameter
// `Problem.EnableInertiaTerms = false`. Have a look at the input file `params.input`
// to see how this is done in this example.
// [[codeblock]]
// We enter the namespace Dumux::Properties in order to import the entire Dumux namespace for general use:
namespace Dumux::Properties {

// declaration of the `ChannelExample` type tag for the single-phase flow problem
namespace TTag {
struct ChannelExample { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // namespace TTag
// [[/codeblock]]

// ### Property specializations
//
// In the following piece of code, mandatory properties for which no meaningful
// default can be set, are specialized for our type tag `ChannelExample`.
// [[codeblock]]
// This sets the grid type used for the simulation. Here, we use a structured 2D grid.
template<class TypeTag>
struct Grid<TypeTag, TTag::ChannelExample> { using type = Dune::YaspGrid<2>; };

// This sets our problem class (see problem.hh) containing initial and boundary conditions.
template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelExample> { using type = Dumux::ChannelExampleProblem<TypeTag> ; };

// This sets the fluid system to be used. Here, we use a liquid with constant properties as fluid phase.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ChannelExample>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};
// [[/codeblock]]

// We also set some properties related to memory management
// throughout the simulation.
// [[details]] caching properties
//
// In Dumux, one has the option to activate/deactivate the grid-wide caching of
// geometries and variables. If active, the CPU time can be significantly reduced
// as less dynamic memory allocation procedures are necessary. Per default, grid-wide
// caching is disabled to ensure minimal memory requirements, however, in this example we
// want to active all available caches, which significantly increases the memory
// demand but makes the simulation faster.
//
// [[codeblock]]
// This enables grid-wide caching of the volume variables.
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };
//This enables grid wide caching for the flux variables.
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };
// This enables grid-wide caching for the finite volume grid geometry
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };
} // end namespace Dumux::Properties
// [[/codeblock]]
// [[/details]]
// [[/content]]
#endif
