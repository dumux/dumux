// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPNCTests
 * \brief The properties of a surfactant model.
 */
#ifndef DUMUX_TEST_2P3C_SURFACTANT_PROPERTIES_HH
#define DUMUX_TEST_2P3C_SURFACTANT_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/2pnc/model.hh>

#include "fluidsystem.hh"
#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

namespace TTag {
struct TestSurfactant
{ using InheritsFrom = std::tuple<TwoPNC>; };

struct TestSurfactantBox
{ using InheritsFrom = std::tuple<TestSurfactant, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TestSurfactant>
{
    using type =  Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TestSurfactant>
{ using type = TestSurfactantProblem<TypeTag>; };

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TestSurfactant>
{ static constexpr bool value = true; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TestSurfactant>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::TestSurfactant<Scalar>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TestSurfactant>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = TestSurfactantSpatialParams<GridGeometry, Scalar>;
};

} // end namespace Dumux::Properties

#endif
