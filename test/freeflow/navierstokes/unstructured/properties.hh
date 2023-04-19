// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties for the channel flow test for the staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_DFG_CHANNEL_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_DFG_CHANNEL_PROPERTIES_HH

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM DFGChannelTestMomentumDiamond
#endif

#ifndef TYPETAG_MASS
#define TYPETAG_MASS DFGChannelTestMassTpfa
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
struct DFGChannelTest {};
struct DFGChannelTestMomentumDiamond { using InheritsFrom = std::tuple<DFGChannelTest, NavierStokesMomentumCVFE, FaceCenteredDiamondModel>; };
struct DFGChannelTestMomentumPQ1Bubble { using InheritsFrom = std::tuple<DFGChannelTest, NavierStokesMomentumCVFE, PQ1BubbleModel>; };
struct DFGChannelTestMassTpfa { using InheritsFrom = std::tuple<DFGChannelTest, NavierStokesMassOneP, CCTpfaModel>; };
struct DFGChannelTestMassBox { using InheritsFrom = std::tuple<DFGChannelTest, NavierStokesMassOneP, BoxModel>; };
struct DFGChannelTestMassDiamond { using InheritsFrom = std::tuple<DFGChannelTest, NavierStokesMassOneP, FaceCenteredDiamondModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DFGChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DFGChannelTest>
{
    using type = Dune::UGGrid<2>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{ using type = DFGChannelTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>> ; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MASS>
{ using type = DFGChannelTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DFGChannelTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DFGChannelTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DFGChannelTest> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DFGChannelTest>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_MOMENTUM, TTag::TYPETAG_MASS>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
