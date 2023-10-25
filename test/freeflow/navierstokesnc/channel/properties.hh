// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesNCTests
 * \brief The properties of the channel flow test for the multi-component staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_CHANNEL_NC_TEST_PROPERTIES_HH
#define DUMUX_CHANNEL_NC_TEST_PROPERTIES_HH

#ifndef ENABLECACHING
#define ENABLECACHING 1
#endif

#ifndef USEMOLES
#define USEMOLES 1
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/fcstaggered.hh>

#include <dumux/freeflow/navierstokes/mass/1pnc/model.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct ChannelNCTest {};
struct ChannelNCTestMomentum { using InheritsFrom = std::tuple<ChannelNCTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
#if !NONISOTHERMAL
struct ChannelNCTestMass { using InheritsFrom = std::tuple<ChannelNCTest, NavierStokesMassOnePNC, CCTpfaModel>; };
#else
struct ChannelNCTestMass { using InheritsFrom = std::tuple<ChannelNCTest, NavierStokesMassOnePNCNI, CCTpfaModel>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelNCTestMomentum>
{ using type = ChannelNCTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelNCTestMass>
{ using type = ChannelNCTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

// Select the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ChannelNCTest>
{
    using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
    static constexpr int phaseIdx = H2OAir::liquidPhaseIdx;
    using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::ChannelNCTest>
{ static constexpr int value = 0; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ChannelNCTest>
{ using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ChannelNCTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ChannelNCTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ChannelNCTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFaceVariablesCache<TypeTag, TTag::ChannelNCTest> { static constexpr bool value = ENABLECACHING; };

// Use mole fraction formulation
template<class TypeTag>
struct UseMoles<TypeTag, TTag::ChannelNCTest> { static constexpr bool value = USEMOLES; };

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::ChannelNCTest>
{
    using Traits = MultiDomainTraits<TTag::ChannelNCTestMomentum, TTag::ChannelNCTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
