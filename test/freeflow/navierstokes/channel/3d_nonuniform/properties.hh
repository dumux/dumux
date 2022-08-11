// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties for the channel flow test for the staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_3D_CHANNEL_PROPERTIES_HH

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM ThreeDChannelTestMomentumDiamond
#endif

#ifndef TYPETAG_MASS
#define TYPETAG_MASS ThreeDChannelTestMassTpfa
#endif

#include <dune/grid/uggrid.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/facecentered/diamond/subcontrolvolumeface.hh>
#include <dumux/discretization/facecentered/diamond/fvgridgeometry.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq1bubble/subcontrolvolumeface.hh>
#include <dumux/discretization/pq1bubble/fvgridgeometry.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/cellcentered/tpfa/subcontrolvolumeface.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/box/subcontrolvolumeface.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>

#include <dumux/freeflow/navierstokes/momentum/diamond/model.hh>
#include <dumux/freeflow/navierstokes/momentum/pq1bubble/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct ThreeDChannelTest {};
struct ThreeDChannelTestMomentumDiamond { using InheritsFrom = std::tuple<ThreeDChannelTest, NavierStokesMomentumDiamond, FaceCenteredDiamondModel>; };
struct ThreeDChannelTestMomentumPQ1Bubble { using InheritsFrom = std::tuple<ThreeDChannelTest, NavierStokesMomentumPQ1Bubble, PQ1BubbleModel>; };
struct ThreeDChannelTestMassTpfa { using InheritsFrom = std::tuple<ThreeDChannelTest, NavierStokesMassOneP, CCTpfaModel>; };
struct ThreeDChannelTestMassBox { using InheritsFrom = std::tuple<ThreeDChannelTest, NavierStokesMassOneP, BoxModel>; };
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
    using type = Dune::UGGrid<3>;
};

// Set the grid type
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::ThreeDChannelTestMomentumDiamond>
{
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();

    // use boundary segment index (works with gmsh files and not with dgf when using ALUGrid)
    struct MyScvfTraits : public FaceCenteredDiamondScvfGeometryTraits<GridView>
    { using BoundaryFlag = BoundarySegmentIndexFlag; };

    struct MyGGTraits : public FaceCenteredDiamondDefaultGridGeometryTraits<GridView>
    { using SubControlVolumeFace = FaceCenteredDiamondSubControlVolumeFace<GridView, MyScvfTraits>; };

    using type = FaceCenteredDiamondFVGridGeometry<GridView, enableCache, MyGGTraits>;
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


template<class TypeTag>
struct GridGeometry<TypeTag, TTag::ThreeDChannelTestMassTpfa>
{
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();

    // use boundary segment index (works with gmsh files and not with dgf when using ALUGrid)
    struct MyScvfTraits : public CCTpfaDefaultScvfGeometryTraits<GridView>
    { using BoundaryFlag = BoundarySegmentIndexFlag; };

    struct MyGGTraits : public CCTpfaDefaultGridGeometryTraits<GridView>
    { using SubControlVolumeFace = CCTpfaSubControlVolumeFace<GridView, MyScvfTraits>; };

    using type = CCTpfaFVGridGeometry<GridView, enableCache, MyGGTraits>;
};

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::ThreeDChannelTestMassBox>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();

    // use boundary segment index (works with gmsh files and not with dgf when using ALUGrid)
    struct MyScvfTraits : public BoxDefaultScvfGeometryTraits<GridView>
    { using BoundaryFlag = BoundarySegmentIndexFlag; };

    struct MyGGTraits : public BoxDefaultGridGeometryTraits<GridView>
    { using SubControlVolumeFace = BoxSubControlVolumeFace<GridView, MyScvfTraits>; };

    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, MyGGTraits>;
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
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
