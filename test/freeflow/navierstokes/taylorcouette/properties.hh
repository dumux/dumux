// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties for the Taylor-Couette benchmark.
 */
#ifndef DUMUX_TAYLOR_COUETTE_PROPERTIES_HH
#define DUMUX_TAYLOR_COUETTE_PROPERTIES_HH

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM TaylorCouetteTestMomentumDiamond
#endif

#ifndef TYPETAG_MASS
#define TYPETAG_MASS TaylorCouetteTestMassTpfa
#endif

#include <dune/grid/uggrid.hh>

#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct TaylorCouetteTest {};
struct TaylorCouetteTestMomentumDiamond { using InheritsFrom = std::tuple<TaylorCouetteTest, NavierStokesMomentumCVFE, FaceCenteredDiamondModel>; };
struct TaylorCouetteTestMomentumPQ1Bubble { using InheritsFrom = std::tuple<TaylorCouetteTest, NavierStokesMomentumCVFE, PQ1BubbleModel>; };
struct TaylorCouetteTestMassTpfa { using InheritsFrom = std::tuple<TaylorCouetteTest, NavierStokesMassOneP, CCTpfaModel>; };
struct TaylorCouetteTestMassBox { using InheritsFrom = std::tuple<TaylorCouetteTest, NavierStokesMassOneP, BoxModel>; };
struct TaylorCouetteTestMassDiamond { using InheritsFrom = std::tuple<TaylorCouetteTest, NavierStokesMassOneP, FaceCenteredDiamondModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TaylorCouetteTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TaylorCouetteTest>
{
    using type = Dune::UGGrid<2>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{ using type = TaylorCouette<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>> ; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MASS>
{ using type = TaylorCouette<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TaylorCouetteTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TaylorCouetteTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TaylorCouetteTest> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TaylorCouetteTest>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_MOMENTUM, TTag::TYPETAG_MASS>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
