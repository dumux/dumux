// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties of the Taylor-Couette test
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_TAYLORCOUETTE_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_TAYLORCOUETTE_PROPERTIES_HH

#include <dune/grid/uggrid.hh>

#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/cctpfa.hh>

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

namespace TTag {
struct TaylorCouetteTest {};
struct TaylorCouetteTestMomentum { using InheritsFrom = std::tuple<TaylorCouetteTest, NavierStokesMomentumCVFE, FaceCenteredDiamondModel>; };
struct TaylorCouetteTestMass { using InheritsFrom = std::tuple<TaylorCouetteTest, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TaylorCouetteTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

template<class TypeTag>
struct Grid<TypeTag, TTag::TaylorCouetteTest>
{ using type = Dune::UGGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TaylorCouetteTestMomentum>
{ using type = TaylorCouetteProblem<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TaylorCouetteTestMass>
{ using type = TaylorCouetteProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::TaylorCouetteTestMomentum> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TaylorCouetteTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TaylorCouetteTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TaylorCouetteTest> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TaylorCouetteTest>
{
    using Traits = MultiDomainTraits<TTag::TaylorCouetteTestMomentum, TTag::TaylorCouetteTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
