// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties of the Darcy-Brinkman test
 */
#ifndef DUMUX_BRINKMAN_TEST_PROPERTIES_HH
#define DUMUX_BRINKMAN_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/components/simpleh2o.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"
#include "spatialparams.hh"

#ifndef NAVIER_STOKES_MODEL
#define NAVIER_STOKES_MODEL NavierStokesMomentum
#endif

#ifndef MOMENTUM_DISCRETIZATION_MODEL
#define MOMENTUM_DISCRETIZATION_MODEL FaceCenteredStaggeredModel
#endif

#ifndef MASS_DISCRETIZATION_MODEL
#define MASS_DISCRETIZATION_MODEL CCTpfaModel
#endif

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>
#endif

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct BrinkmanTest {};
struct BrinkmanTestMomentum { using InheritsFrom = std::tuple<BrinkmanTest, NAVIER_STOKES_MODEL, MOMENTUM_DISCRETIZATION_MODEL>; };
struct BrinkmanTestMass { using InheritsFrom = std::tuple<BrinkmanTest, NavierStokesMassOneP, MASS_DISCRETIZATION_MODEL>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::BrinkmanTestMomentum>
{ using type = BrinkmanProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::BrinkmanTestMass>
{ using type = BrinkmanProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::BrinkmanTest>
{ using type = BrinkmanTestSpatialParams<GetPropType<TypeTag, GridGeometry>, GetPropType<TypeTag, Scalar>>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::BrinkmanTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::BrinkmanTest>
{ using type = GRIDTYPE; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::BrinkmanTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::BrinkmanTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::BrinkmanTest>
{ static constexpr bool value = true; };

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::BrinkmanTest>
{
    using Traits = MultiDomainTraits<TTag::BrinkmanTestMomentum, TTag::BrinkmanTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
