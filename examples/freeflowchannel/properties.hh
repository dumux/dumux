// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_EXAMPLES_FREEFLOW_CHANNEL_PROPERTIES_HH
#define DUMUX_EXAMPLES_FREEFLOW_CHANNEL_PROPERTIES_HH

// ## Compile-time settings (`properties.hh`)
//
// In this file, the type tag used for this simulation (`TTag::ChannelExample`) is defined.
// As this is a coupled problem we also define type tags for the two subproblems
// `TTag::ChannelExampleMass` and `TTag::ChannelExampleMomentum`. We then specialize properties
// (compile time options) to the needs of the desired setup for the respective type tags.
//
// [[content]]
//
// ### Includes
// [[details]] includes
//
// The `NavierStokesMomentum` and `NavierStokesMass` type tags specialize most of the properties
// required for Navier-Stokes single-phase flow simulations in DuMu<sup>x</sup>. We will use this in
// the following to inherit the respective properties and subsequently specialize those properties
// for our type tags, which we want to modify or for which no meaningful default can be set.
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>

// We want to use `YaspGrid`, an implementation of the dune grid interface for structured grids:
#include <dune/grid/yaspgrid.hh>
// In this example, we want to discretize the momentum and mass balances with the staggered-grid and
// cell-centered finite volume discretization schemes, respectively:
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
// The fluid properties are specified in the following headers (we use a liquid with constant properties as the fluid phase):
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

// We include the problem header used for this simulation.
#include "problem.hh"
// [[/details]]
//
// ### Type tag definition
//
// We define a type tag for our simulation with the name `ChannelExample` as well as type tags for
// the momentum and mass subproblems. Shared properties can be specialized to `ChannelExample` and
// will be inherited by the subproblems' type tags.
// These inherit the properties specialized for the physical models `NavierStokesMomentum` and
// `NavierStokesMassOneP` as well as the discretization methods `FaceCenteredStaggeredModel`
// and `CCTpfaModel` respectively.
// This way, most of the properties required for Navier-Stokes single-phase flow simulations
// using the staggered-grid scheme are conveniently specialized for our new type tag.
// However, some properties depend on user choices and no meaningful default value can be set.
// Those properties will be addressed later in this file.
// Please note that, in this example, we actually want to solve the Stokes instead of the
// Navier-Stokes equations. This can be achieved at runtime by setting the parameter
// `Problem.EnableInertiaTerms = false`. Have a look at the input file `params.input`
// to see how this is done in this example.
// [[codeblock]]
// We enter the namespace `Dumux::Properties` in order to import the entire `Dumux` namespace for general use:
namespace Dumux::Properties {

// declaration of the `ChannelExample` type tag for the single-phase flow problem
namespace TTag {
struct ChannelExample {};
struct ChannelExampleMomentum { using InheritsFrom = std::tuple<ChannelExample, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct ChannelExampleMass { using InheritsFrom = std::tuple<ChannelExample, NavierStokesMassOneP, CCTpfaModel>; };
} // namespace TTag
// [[/codeblock]]

// ### Property specializations
//
// In the following piece of code, mandatory properties for which no meaningful default can be set,
// are specialized for our type tag `ChannelExample` or the appropriate type tag of a subproblem.
// [[codeblock]]
// This sets the grid type used for the simulation. Here, we use a structured 2D grid.
template<class TypeTag>
struct Grid<TypeTag, TTag::ChannelExample> { using type = Dune::YaspGrid<2>; };

// This sets our problem type (see `problem.hh`) containing the initial and boundary conditions.
template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelExampleMomentum>
{ using type = Dumux::ChannelExampleProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>> ; };
template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelExampleMass>
{ using type = Dumux::ChannelExampleProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>> ; };

// This sets the fluid system type to be used. Here, we use a liquid with constant properties as the fluid phase.
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
// In DuMu<sup>x</sup>, one has the option to activate/deactivate the grid-wide caching of
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
// [[/codeblock]]
// [[/details]]

// Finally we define the coupling manager to couple the momentum and mass subproblems
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::ChannelExample>
{
    using Traits = MultiDomainTraits<TTag::ChannelExampleMomentum, TTag::ChannelExampleMass>;
    using type = FreeFlowCouplingManager<Traits>;
};
} // end namespace Dumux::Properties
// [[/content]]
#endif
