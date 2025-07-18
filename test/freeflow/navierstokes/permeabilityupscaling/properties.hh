// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties of the pore flow test for permeability upscaling
 */
#ifndef DUMUX_UPSCALING_TEST_PROPERTIES_HH
#define DUMUX_UPSCALING_TEST_PROPERTIES_HH

#ifndef ENABLECACHING
#define ENABLECACHING 1
#endif

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM PoreFlowTestMomentumStaggered
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/subgrid/subgrid.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"
namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct PoreFlowTest {};
struct PoreFlowTestMomentumStaggered { using InheritsFrom = std::tuple<PoreFlowTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct PoreFlowTestMomentumDiamond { using InheritsFrom = std::tuple<PoreFlowTest, NavierStokesMomentumCVFE, FaceCenteredDiamondModel>; };
struct PoreFlowTestMass { using InheritsFrom = std::tuple<PoreFlowTest, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{ using type = PoreFlowTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::PoreFlowTestMass>
{ using type = PoreFlowTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::PoreFlowTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::PoreFlowTest>
{
    static constexpr int dim = 3;
    using HostGrid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, dim> >;
    using type = Dune::SubGrid<HostGrid::dimension, HostGrid>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::PoreFlowTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::PoreFlowTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::PoreFlowTest> { static constexpr bool value = ENABLECACHING; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::PoreFlowTest>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_MOMENTUM, TTag::PoreFlowTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
