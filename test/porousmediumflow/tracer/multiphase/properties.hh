// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup TracerTests
 * \brief The test properties
 */

#ifndef DUMUX_TRACER_MULTIPHASE_TEST_PROPERTIES_HH
#define DUMUX_TRACER_MULTIPHASE_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/porousmediumflow/tracer/model.hh>

#include "fluidsystem.hh"
#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct TracerTest { using InheritsFrom = std::tuple<Tracer>; };
struct TracerTestTpfa { using InheritsFrom = std::tuple<TracerTest, CCTpfaModel>; };
struct TracerTestMpfa { using InheritsFrom = std::tuple<TracerTest, CCMpfaModel>; };
struct TracerTestBox { using InheritsFrom = std::tuple<TracerTest, BoxModel>; };
} // end namespace TTag

// enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TracerTest> { static constexpr bool value = true; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TracerTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TracerTest> { using type = TracerTest<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TracerTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TracerTestSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::TracerTest> { static constexpr bool value = false; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TracerTest> { using type = FluidSystems::TracerTest<GetPropType<TypeTag, Properties::Scalar>>; };

} // end namespace Dumux::Properties

#endif
