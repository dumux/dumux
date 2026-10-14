// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Properties of the lock-exchange test with the Boussinesq approximation
 */
#ifndef DUMUX_LOCKEXCHANGE_TEST_PROBLEM_PROPERTIES_HH
#define DUMUX_LOCKEXCHANGE_TEST_PROBLEM_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box.hh>
#include <dumux/flux/cvfe/boussinesqdarcyslaw.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/liquidphase2c.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>

#include "problem.hh"
#include "../../spatialparams.hh"

namespace Dumux::Properties {

namespace TTag {
struct LockExchangeBoussinesq { using InheritsFrom = std::tuple<OnePNC, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::LockExchangeBoussinesq> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::LockExchangeBoussinesq> { using type = LockExchangeTestProblem<TypeTag>; };

// components with equal molar masses give a constant (reference) density
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::LockExchangeBoussinesq>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::LiquidPhaseTwoC<Scalar, Components::Constant<1, Scalar>, Components::Constant<2, Scalar>>;
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::LockExchangeBoussinesq>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNCTestSpatialParams<GridGeometry, Scalar>;
};

template<class TypeTag>
struct UseMoles<TypeTag, TTag::LockExchangeBoussinesq> { static constexpr bool value = false; };

template<class TypeTag>
struct AdvectionType<TypeTag, TTag::LockExchangeBoussinesq>
{
    using type = CVFEBoussinesqDarcysLaw<GetPropType<TypeTag, Properties::Scalar>,
                                         GetPropType<TypeTag, Properties::GridGeometry>>;
};

} // end namespace Dumux::Properties

#endif
