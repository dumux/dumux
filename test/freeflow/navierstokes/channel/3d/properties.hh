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
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROPERTIES_HH

#ifndef ENABLECACHING
#define ENABLECACHING 1
#endif

#ifndef GRID_DIM
#define GRID_DIM 3
#endif

#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#endif

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "spatialparams.hh"
#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct ThreeDChannelTest {};
struct ThreeDChannelTestMomentum { using InheritsFrom = std::tuple<ThreeDChannelTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct ThreeDChannelTestMass { using InheritsFrom = std::tuple<ThreeDChannelTest, NavierStokesMassOneP, CCTpfaModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ThreeDChannelTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ThreeDChannelTest>
{
    static constexpr int dim = GRID_DIM;

    using HostGrid = Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, dim> >;

#if HAVE_DUNE_SUBGRID
    using type = Dune::SubGrid<HostGrid::dimension, HostGrid>;
#else
    using type = HostGrid;
#endif
};

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::ThreeDChannelTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = Channel3DSpatialParams<GridGeometry, Scalar>;
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ThreeDChannelTestMomentum>
{ using type = ThreeDChannelTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::ThreeDChannelTestMass>
{ using type = ThreeDChannelTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ThreeDChannelTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ThreeDChannelTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ThreeDChannelTest> { static constexpr bool value = ENABLECACHING; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::ThreeDChannelTest>
{
    using Traits = MultiDomainTraits<TTag::ThreeDChannelTestMomentum, TTag::ThreeDChannelTestMass>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
