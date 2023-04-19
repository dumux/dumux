// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Definition of a problem for the 1pnc problem:
 * Component transport of nitrogen dissolved in the water phase.
 */

#ifndef DUMUX_1P2C_TEST_PROBLEM_PROPERTIES_HH
#define DUMUX_1P2C_TEST_PROBLEM_PROPERTIES_HH

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/math.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>


#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include "problem.hh"
#include "../spatialparams.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct OnePTwoCTest { using InheritsFrom = std::tuple<OnePNC>; };
struct OnePTwoCTestBox { using InheritsFrom = std::tuple<OnePTwoCTest, BoxModel>; };
struct OnePTwoCTestCCTpfa { using InheritsFrom = std::tuple<OnePTwoCTest, CCTpfaModel>; };
struct OnePTwoCTestCCMpfa { using InheritsFrom = std::tuple<OnePTwoCTest, CCMpfaModel>; };
} // end namespace TTag

// Set the grid type
#if HAVE_DUNE_UGGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePTwoCTest> { using type = Dune::UGGrid<2>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePTwoCTest> { using type = Dune::YaspGrid<2>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePTwoCTest> { using type = OnePTwoCTestProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePTwoCTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using H2ON2 = FluidSystems::H2ON2<Scalar, FluidSystems::H2ON2DefaultPolicy</*simplified=*/true>>;
    using type = FluidSystems::OnePAdapter<H2ON2, H2ON2::liquidPhaseIdx>;
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePTwoCTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNCTestSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::OnePTwoCTest> { static constexpr bool value = true; };
} // end namespace Dumux::Properties

#endif
