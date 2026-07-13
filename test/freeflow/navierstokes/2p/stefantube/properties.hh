// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Properties for the Stefan-tube slit-capillary evaporation benchmark.
 */
#ifndef DUMUX_TEST_FREEFLOW_TWOP_STEFAN_TUBE_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_TWOP_STEFAN_TUBE_PROPERTIES_HH

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM StefanTubeMomentumPQ1Bubble
#endif

#ifndef TYPETAG_MASS
#define TYPETAG_MASS StefanTubeMassBox
#endif

#ifndef STEFANTUBE_DIM
#define STEFANTUBE_DIM 2
#endif

#include <dune/alugrid/grid.hh>

#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/box.hh>

#include <dumux/freeflow/navierstokes/mass/2pvapor/model.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/2p/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"
#include "../couplingmanager.hh"

namespace Dumux::Properties {

namespace TTag {
struct StefanTube {};
struct StefanTubeMomentumPQ1Bubble { using InheritsFrom = std::tuple<StefanTube, NavierStokesMomentumTwoPCVFE, PQ1BubbleModel>; };
struct StefanTubeMassBox { using InheritsFrom = std::tuple<StefanTube, NavierStokesMassTwoPVapor, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::StefanTube>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar>>;
};

template<class TypeTag>
struct FluidState<TypeTag, TTag::StefanTube>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

template<class TypeTag>
struct Grid<TypeTag, TTag::StefanTube>
{ using type = Dune::ALUGrid<STEFANTUBE_DIM, STEFANTUBE_DIM, Dune::simplex, Dune::conforming>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{ using type = StefanTubeProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MASS>
{ using type = StefanTubeProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::StefanTube> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::StefanTube> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::StefanTube> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::StefanTube>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_MOMENTUM, TTag::TYPETAG_MASS>;
    using type = CVFEFreeFlowCouplingManagerTwoP<Traits>;
};

} // end namespace Dumux::Properties

#endif
