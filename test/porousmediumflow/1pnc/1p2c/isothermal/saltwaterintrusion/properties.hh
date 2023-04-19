// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Definition of a problem involving salt
 *        water intrusion into a fresh water aquifer.
 */

#ifndef DUMUX_SALTWATERINTRUSION_TEST_PROBLEM_PROPERTIES_HH
#define DUMUX_SALTWATERINTRUSION_TEST_PROBLEM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>


#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include <dumux/material/fluidsystems/brine.hh>

#include "problem.hh"
#include "../../spatialparams.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct SaltWaterIntrusionTest { using InheritsFrom = std::tuple<OnePNC, BoxModel>; };
} // end namespace TTag

// Use a structured yasp grid
template<class TypeTag>
struct Grid<TypeTag, TTag::SaltWaterIntrusionTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::SaltWaterIntrusionTest> { using type = SaltWaterIntrusionTestProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::SaltWaterIntrusionTest>
{ using type = FluidSystems::Brine< GetPropType<TypeTag, Properties::Scalar> >; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::SaltWaterIntrusionTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNCTestSpatialParams<GridGeometry, Scalar>;
};

// Use mass fractions to set salinity conveniently
template<class TypeTag>
struct UseMoles<TypeTag, TTag::SaltWaterIntrusionTest> { static constexpr bool value = false; };

} // end namespace Dumux::Properties

#endif
