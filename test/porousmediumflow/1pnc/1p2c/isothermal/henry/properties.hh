// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief The Henry problem benchmark of Fahs et al. (2016, WRR,
 *        doi:10.1002/2016WR019288): property definitions.
 */

#ifndef DUMUX_HENRY_FAHS_TEST_PROBLEM_PROPERTIES_HH
#define DUMUX_HENRY_FAHS_TEST_PROBLEM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "problem.hh"
#include "fluidsystem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct HenryFahsTest { using InheritsFrom = std::tuple<OnePNC, BoxModel>; };
} // end namespace TTag

// Use a structured yasp grid
template<class TypeTag>
struct Grid<TypeTag, TTag::HenryFahsTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::HenryFahsTest> { using type = HenryFahsTestProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::HenryFahsTest>
{ using type = FluidSystems::HenryFahsFluid<GetPropType<TypeTag, Properties::Scalar>>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::HenryFahsTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = HenryFahsSpatialParams<GridGeometry, Scalar>;
};

// Use mass fractions to set salinity conveniently
template<class TypeTag>
struct UseMoles<TypeTag, TTag::HenryFahsTest> { static constexpr bool value = false; };

} // end namespace Dumux::Properties

#endif
