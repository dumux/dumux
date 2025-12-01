// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties for the channel flow test for the staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROPERTIES_HH

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM ThreeDChannelTestMomentumPQ1Bubble
#endif

#ifndef TYPETAG_MASS
#define TYPETAG_MASS ThreeDChannelTestMassBox
#endif

#include <dune/alugrid/grid.hh>

#include <dumux/common/boundaryflag.hh>

#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq1bubble/subcontrolvolumeface.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/box/subcontrolvolumeface.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>

#include <dumux/freeflow/navierstokes/mass/2p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/2p/cvfe/model.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/2pimmiscible.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"
#include "couplingmanager.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct ThreeDChannelTest {};
struct ThreeDChannelTestMomentumPQ1Bubble { using InheritsFrom = std::tuple<ThreeDChannelTest, NavierStokesMomentumTwoPCVFE, PQ1BubbleModel>; };
struct ThreeDChannelTestMassBox { using InheritsFrom = std::tuple<ThreeDChannelTest, NavierStokesMassTwoP, BoxModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ThreeDChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar>>;
};

template<class TypeTag>
struct FluidState<TypeTag, TTag::ThreeDChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ThreeDChannelTest>
{
    using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>;
};

// Set the grid type
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::ThreeDChannelTestMomentumPQ1Bubble>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();

    // use boundary segment index (works with gmsh files and not with dgf when using ALUGrid)
    struct MyScvfTraits : public PQ1BubbleDefaultScvfGeometryTraits<GridView>
    { using BoundaryFlag = BoundarySegmentIndexFlag; };

    struct MyGGTraits : public PQ1BubbleDefaultGridGeometryTraits<GridView>
    { using SubControlVolumeFace = PQ1BubbleSubControlVolumeFace<GridView, MyScvfTraits>; };

    using type = PQ1BubbleFVGridGeometry<Scalar, GridView, enableCache, MyGGTraits>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{ using type = ThreeDChannelTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MASS>
{ using type = ThreeDChannelTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ThreeDChannelTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ThreeDChannelTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ThreeDChannelTest> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::ThreeDChannelTest>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_MOMENTUM, TTag::TYPETAG_MASS>;
    using type = CVFEFreeFlowCouplingManagerTwoP<Traits>;
};

} // end namespace Dumux::Properties

#endif
