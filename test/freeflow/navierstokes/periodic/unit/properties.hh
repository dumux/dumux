// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Stationary test for the staggered grid Navier-Stokes model with periodic BC
 */

#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_PERIODIC_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_PERIODIC_PROPERTIES_HH

#ifndef MOMENTUM_TYPETAG
#define MOMENTUM_TYPETAG PeriodicTestMomentumStaggered
#endif

#ifndef MASS_TYPETAG
#define MASS_TYPETAG PeriodicTestMassTpfa
#endif

#include <dune/grid/spgrid.hh>
#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>

#include "problem.hh"

#ifndef USESUBGRID
#define USESUBGRID 0
#endif

#ifndef USEALUGRID
#define USEALUGRID 0
#endif

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PeriodicTest {};
struct PeriodicTestMomentumStaggered { using InheritsFrom = std::tuple<PeriodicTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct PeriodicTestMomentumDiamond { using InheritsFrom = std::tuple<PeriodicTest, NavierStokesMomentumCVFE, FaceCenteredDiamondModel>; };
struct PeriodicTestMomentumPQ1Bubble { using InheritsFrom = std::tuple<PeriodicTest, NavierStokesMomentumCVFE, PQ1BubbleModel>; };
struct PeriodicTestMassTpfa { using InheritsFrom = std::tuple<PeriodicTest, NavierStokesMassOneP, CCTpfaModel>; };
struct PeriodicTestMassBox { using InheritsFrom = std::tuple<PeriodicTest, NavierStokesMassOneP, BoxModel>; };
struct PeriodicTestMassDiamond { using InheritsFrom = std::tuple<PeriodicTest, NavierStokesMassOneP, FaceCenteredDiamondModel>; };
} // end namespace TTag


// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PeriodicTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PeriodicTest> {
#if HAVE_DUNE_ALUGRID && USEALUGRID
    using HostGrid = Dune::ALUGrid<2,2,Dune::simplex,Dune::nonconforming>;
#else
    using HostGrid = Dune::SPGrid<double, 2>;
#endif
#if HAVE_DUNE_SUBGRID && USESUBGRID
    using type = Dune::SubGrid<2, HostGrid>;
#else
    using type = HostGrid;
#endif
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::MOMENTUM_TYPETAG>
{ using type = PeriodicTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::MASS_TYPETAG>
{ using type = PeriodicTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PeriodicTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::PeriodicTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PeriodicTest> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PeriodicTest>
{
    using Traits = MultiDomainTraits<TTag::MOMENTUM_TYPETAG, TTag::MASS_TYPETAG>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
