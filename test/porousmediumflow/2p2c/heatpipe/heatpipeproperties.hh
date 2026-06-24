// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux-Lecture contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_HEATPIPE_PROPERTIES_HH
#define DUMUX_HEATPIPE_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>

#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/porousmediumflow/2p2c/model.hh>

#include "heatpipespatialparams.hh"
#include "heatpipeproblem.hh"

namespace Dumux::Properties
{
// Create new type tags
namespace TTag {
struct HeatPipeTypeTag { using InheritsFrom = std::tuple<TwoPTwoCNI, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::HeatPipeTypeTag> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::HeatPipeTypeTag> { using type = HeatPipeProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::HeatPipeTypeTag> { using type = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>; };

// pn-sw formulation
template<class TypeTag>
struct Formulation<TypeTag, TTag::HeatPipeTypeTag> { static constexpr auto value = TwoPFormulation::p1s0; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::HeatPipeTypeTag>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = HeatPipeSpatialParams<FVGridGeometry, Scalar>;
};

} // namespace Dumux::Properties

#endif
