// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief The properties of the non-isothermal gas injection problem where a gas (e.g. air) is injected into a fully water saturated medium.
 */

#ifndef DUMUX_WATER_AIR_PROPERTIES_HH
#define DUMUX_WATER_AIR_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/porousmediumflow/2p2c/model.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct WaterAir { using InheritsFrom = std::tuple<TwoPTwoCNI>; };
struct WaterAirBox { using InheritsFrom = std::tuple<WaterAir, BoxModel>; };
struct WaterAirCCTpfa { using InheritsFrom = std::tuple<WaterAir, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::WaterAir> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::WaterAir> { using type = WaterAirProblem<TypeTag>; };

// Set the wetting phase
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::WaterAir> { using type = FluidSystems::H2ON2<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::WaterAir>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = WaterAirSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };

// Enable caching
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::WaterAir> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
