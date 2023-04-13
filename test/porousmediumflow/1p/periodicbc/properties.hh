// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The properties for the incompressible test
 */

#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_PROPERTIES_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_PROPERTIES_HH

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#endif

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"
#include "spatialparams.hh"

#ifndef FVGEOMCACHING
#define FVGEOMCACHING 0
#endif

namespace Dumux::Properties {
namespace TTag {
struct OnePIncompressible { using InheritsFrom = std::tuple<OneP>; };
struct OnePIncompressibleTpfa { using InheritsFrom = std::tuple<OnePIncompressible, CCTpfaModel>; };
struct OnePIncompressibleMpfa { using InheritsFrom = std::tuple<OnePIncompressible, CCMpfaModel>; };
struct OnePIncompressibleBox { using InheritsFrom = std::tuple<OnePIncompressible, BoxModel>; };
} // end namespace TTag

// Set the grid type
#if HAVE_DUNE_SPGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePIncompressible> { using type = Dune::SPGrid<double, 2>; };
#endif

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePIncompressible> { using type = OnePTestProblem<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePIncompressible> { using type = OnePTestSpatialParams<TypeTag>; };

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePIncompressible> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePIncompressible>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = FVGEOMCACHING; };
} // end namespace Dumux::Properties
#endif
