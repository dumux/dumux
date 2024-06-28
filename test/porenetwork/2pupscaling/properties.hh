// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_PNM_UPSCALING_PROPERTIES_HH
#define DUMUX_PNM_UPSCALING_PROPERTIES_HH
/*!
 * \file
 *
 * \brief Test for the pore network model
 */
 #include <config.h>

 #include <ctime>
 #include <iostream>

 #include <dune/common/parallel/mpihelper.hh>
 #include <dune/common/timer.hh>
 #include <dune/grid/io/file/dgfparser/dgfexception.hh>
 #include <dune/grid/io/file/vtk.hh>
 #include <dune/grid/io/file/vtk/vtksequencewriter.hh>

 #include <dumux/common/initialize.hh>
 #include <dumux/common/properties.hh>
 #include <dumux/common/parameters.hh>
 #include <dumux/common/dumuxmessage.hh>

#include <dumux/common/properties/model.hh>
#include <dumux/common/properties/grid.hh>
#include <dumux/discretization/porenetwork/gridgeometry.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/thresholdcapillarypressures.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/pore/2p/localrulesforplatonicbody.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility1p.hh>
#include <dumux/material/fluidmatrixinteractions/porenetwork/throat/transmissibility2p.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

#include <dumux/porenetwork/2p/static/staticdrainge.hh>
#include <dumux/io/gnuplotinterface.hh>

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
#include <dumux/material/fluidsystems/1pgas.hh>

#include "helper.hh"
#include "problem_static.hh"
#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PNMTWOPStatic { using InheritsFrom = std::tuple<GridProperties, ModelProperties>; };
struct PNMUpscalingCreepingFlowLiquid { using InheritsFrom = std::tuple<PNMOneP>; };
struct PNMUpscalingCreepingFlowGas { using InheritsFrom = std::tuple<PNMUpscalingCreepingFlowLiquid>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMTWOPStatic> { using type = Dune::FoamGrid<1, 3>; };

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::PNMTWOPStatic>
{
private:
    static constexpr bool enableCache = false;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = Dumux::PoreNetwork::GridGeometry<Scalar, GridView, enableCache>;
};

// The flux variables cache
template<class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::PNMTWOPStatic>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::SimpleFluxVariablesCache<Scalar>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMTWOPStatic> { using type = Dumux::PoreNetwork::DrainageProblemStatic<TypeTag>; };


// ### Property specializations
//
// In the following piece of code, mandatory `properties` for which no meaningful
// default can be set, are specialized for our type tag `PNMUpscaling`.
// [[codeblock]]
// We use `dune-foamgrid`, which is especially tailored for 1D networks.
template<class TypeTag>
struct Grid<TypeTag, TTag::PNMUpscalingCreepingFlowLiquid>
{ using type = Dune::FoamGrid<1, 3>; };

// The problem class specifying initial and boundary conditions:
template<class TypeTag>
struct Problem<TypeTag, TTag::PNMUpscalingCreepingFlowLiquid>
{ using type = UpscalingProblem<TypeTag>; };

//! The spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::PNMUpscalingCreepingFlowLiquid>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = PoreNetwork::UpscalingSpatialParams<GridGeometry, Scalar>;
};

//! The advection type for creeping flow
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::PNMUpscalingCreepingFlowLiquid>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using TransmissibilityLaw = PoreNetwork::SinglePhaseTransmissibility<Scalar>;
public:
    using type = PoreNetwork::CreepingFlow<Scalar, TransmissibilityLaw>;
};

// We use a single liquid phase consisting of a component with constant fluid properties.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMUpscalingCreepingFlowLiquid>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// We use a single liquid phase consisting of a component with constant fluid properties.
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PNMUpscalingCreepingFlowGas>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePGas<Scalar, Components::Constant<1, Scalar> >;
};
// [[/codeblock]]

// Moreover, here we use a local residual specialized for incompressible flow
// that contains functionality related to analytic differentiation.
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::PNMUpscalingCreepingFlowLiquid>
{ using type = OnePIncompressibleLocalResidual<TypeTag>; };


} // end namespace Dumux::Properties
#endif
