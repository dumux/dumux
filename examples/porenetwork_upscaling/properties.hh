// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_PROPERTIES_HH
#define DUMUX_PNM_ONEP_PERMEABILITY_UPSCALING_PROPERTIES_HH

// ## Compile-time settings (`properties.hh`)
// This file defines the `TypeTag` used for the simulation in this example, for
// which we specialize a number of compile-time `properties`.
// [[content]]
// ### Includes
// [[details]] includes
#include <dune/foamgrid/foamgrid.hh> // for `Dune::FoamGrid`

// The `OneP` type tag specializes most of the `properties` required for single-
// phase flow simulations in DuMu<sup>x</sup>. We will use this in the following to inherit the
// respective properties, and subsequently specialize those properties for our
// type tag, which we want to modify or for which no meaningful default can be set.
#include <dumux/porenetwork/1p/model.hh>// for `TTag::PNMOneP`

// The class that contains a collection of single-phase flow throat transmissibilities
// among them the transmisibility model to be used can be specified in AdvectionType class
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>

// The class that provides specializations for both creeping and non-creeping advection types.
#include <dumux/flux/porenetwork/advection.hh>

// The local residual for incompressible flow is included.
// The one-phase flow model (included above) uses a default implementation of the
// local residual for single-phase flow. However, in this example we are using an
// incompressible fluid phase. Therefore, we are including the specialized local
// residual which contains functionality to analytically compute the entries of
// the Jacobian matrix. We will use this in the main file.
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

// We will use a single liquid phase consisting of a component with constant fluid properties.
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

// The classes that define the problem and parameters used in this simulation
#include "problem.hh"
#include "spatialparams.hh"
// [[/details]]
//
// ### `TypeTag` definition
// Two `TypeTag` for our simulation are defined, one for creeping flow and another for non-creeping flow,
// which inherit properties from the single-phase pore network model. The non-creeping flow inherits
// all properties from the creeping flow simulation but overrides the `AdvectionType` property.
namespace Dumux::Properties {
namespace TTag {
struct PNMUpscalingCreepingFlow { using InheritsFrom = std::tuple<PNMOneP>; };
struct PNMUpscalingNonCreepingFlow { using InheritsFrom = std::tuple<PNMUpscalingCreepingFlow>; };
}

// ### Property specializations
//
// In the following piece of code, mandatory `properties` for which no meaningful
// default can be set, are specialized for our type tag `PNMUpscaling`.
// [[codeblock]]
// We use `dune-foamgrid`, which is especially tailored for 1D networks.
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMUpscalingCreepingFlow>
{ using type = Dune::FoamGrid<1, 3>; };

// The problem class specifying initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMUpscalingCreepingFlow>
{ using type = UpscalingProblem<TypeTag>; };

//! The spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMUpscalingCreepingFlow>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::UpscalingSpatialParams<GridGeometry, Scalar>;
};

//! The advection type for creeping flow
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMUpscalingCreepingFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::TransmissibilityPatzekSilin<Scalar, true/*considerPoreBodyResistance*/>;
public:
    using type = PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

//! The advection type for non-creeping flow (includes model for inertia effects)
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMUpscalingNonCreepingFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::TransmissibilityPatzekSilin<Scalar, true/*considerPoreBodyResistance*/>;
public:
    using type = PoreNetwork::NonCreepingFlow<Scalar, TransmissibilityLaw>;
};

// We use a single liquid phase consisting of a component with constant fluid properties.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMUpscalingCreepingFlow>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};
// [[/codeblock]]

// Moreover, here we use a local residual specialized for incompressible flow
// that contains functionality related to analytic differentiation.
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PNMUpscalingCreepingFlow>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };

} // end namespace Dumux::Properties
// [[/content]]
#endif
