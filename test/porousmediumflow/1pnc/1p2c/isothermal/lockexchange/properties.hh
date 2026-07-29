// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Properties for the lock-exchange test: a "full" (variable-density) type tag
 *        and a Boussinesq-approximated one, sharing the same problem/grid/BCs/IC, so
 *        the two can be regression-tested against each other.
 */

#ifndef DUMUX_LOCKEXCHANGE_TEST_PROBLEM_PROPERTIES_HH
#define DUMUX_LOCKEXCHANGE_TEST_PROBLEM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/flux/cvfe/boussinesqdarcyslaw.hh>

#include "problem.hh"
#include "fluidsystem.hh"
#include "../../spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct LockExchangeFullTest { using InheritsFrom = std::tuple<OnePNC, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::LockExchangeFullTest> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::LockExchangeFullTest> { using type = LockExchangeTestProblem<TypeTag>; };

// the "full" model: density varies with the tracer mass fraction everywhere
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LockExchangeFullTest>
{ using type = FluidSystems::LockExchangeFluid<GetPropType<TypeTag, Properties::Scalar>, true>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LockExchangeFullTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNCTestSpatialParams<GridGeometry, Scalar>;
};

// Use mass fractions to set the tracer concentration conveniently
template<class TypeTag>
struct UseMoles<TypeTag, TTag::LockExchangeFullTest> { static constexpr bool value = false; };

// Boussinesq-approximated variant: same grid, BCs, IC and spatial params, but the
// fluid system's density is the constant reference (tracer-independent), and
// AdvectionType supplies the tracer-driven buoyancy separately via
// BoussinesqCVFEDarcyLaw + Problem::buoyantDensity().
namespace TTag {
struct LockExchangeBoussinesqTest { using InheritsFrom = std::tuple<LockExchangeFullTest>; };
} // end namespace TTag

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LockExchangeBoussinesqTest>
{ using type = FluidSystems::LockExchangeFluid<GetPropType<TypeTag, Properties::Scalar>, false>; };

template<class TypeTag>
struct AdvectionType<TypeTag, TTag::LockExchangeBoussinesqTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using type = BoussinesqCVFEDarcyLaw<Scalar, GridGeometry>;
};

} // end namespace Dumux::Properties

#endif
