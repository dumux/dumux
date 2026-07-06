// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Properties for the P2-P1-P2-P2 Taylor-Hood Cahn-Hilliard / Navier-Stokes rising-bubble test.
 *
 * Two subdomains: (0) co-located momentum + Cahn-Hilliard on PQ2 (Taylor-Hood P2 velocity, P2 c,
 * P2 mu; numEq = dim+2), (1) pressure/continuity on Box P1 (reused single-phase mass model). The
 * density contrast enters only the momentum buoyancy (rho(c) local in the momentum volvars), so the
 * pressure/continuity is the ordinary incompressible one and the standard FreeFlowCouplingManager
 * (momentum <-> pressure) applies unchanged.
 */
#ifndef DUMUX_TEST_FF_NAVIERSTOKES_2P_P2_PROPERTIES_HH
#define DUMUX_TEST_FF_NAVIERSTOKES_2P_P2_PROPERTIES_HH

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM RisingBubbleP2Momentum
#endif
#ifndef TYPETAG_MASS
#define TYPETAG_MASS RisingBubbleP2Mass
#endif

#include <dune/alugrid/grid.hh>

#include <dumux/common/boundaryflag.hh>

#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/pq2/subcontrolvolumeface.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/gridvariables.hh>

#include <dumux/freeflow/navierstokes/momentum/chns/cvfe/model.hh>
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
struct RisingBubbleP2 {};
struct RisingBubbleP2Momentum { using InheritsFrom = std::tuple<RisingBubbleP2, NavierStokesMomentumCHNSCVFE, PQ2HybridModel>; };
struct RisingBubbleP2Mass { using InheritsFrom = std::tuple<RisingBubbleP2, NavierStokesMassOneP, BoxModel>; };
} // end namespace TTag

// fluid system (single-phase constant; the two-phase material is provided by the problem)
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::RisingBubbleP2>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<0, Scalar>>;
};

// grid: 2D simplex, conforming (P2 space stays conforming; matches the rect.msh gmsh mesh)
template<class TypeTag>
struct Grid<TypeTag, TTag::RisingBubbleP2>
{ using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };

// PQ2 (Taylor-Hood) momentum grid geometry with gmsh boundary-segment flags
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::RisingBubbleP2Momentum>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();

    struct MyScvfTraits : public PQ2DefaultScvfGeometryTraits<GridView>
    { using BoundaryFlag = BoundarySegmentIndexFlag; };

    struct MyGGTraits : public PQ2DefaultGridGeometryTraits<GridView>
    { using SubControlVolumeFace = PQ2SubControlVolumeFace<GridView, MyScvfTraits>; };

    using type = PQ2FVGridGeometry<Scalar, GridView, enableCache, MyGGTraits>;
};

// hybrid CVFE grid variables (P2 = vertex CV + FE edge dofs) for the momentum+CH subdomain
template<class TypeTag>
struct GridVariables<TypeTag, TTag::RisingBubbleP2Momentum>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

// experimental (new-interface) grid variables for the P1 Box pressure/mass subdomain
template<class TypeTag>
struct GridVariables<TypeTag, TTag::RisingBubbleP2Mass>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::CVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

// problems: co-located momentum+CH (CVFE momentum problem base) and pressure (CVFE mass problem base)
template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{ using type = RisingBubbleP2Problem<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MASS>
{ using type = RisingBubbleP2Problem<TypeTag, Dumux::CVFENavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::RisingBubbleP2> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::RisingBubbleP2> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::RisingBubbleP2> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::RisingBubbleP2>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_MOMENTUM, TTag::TYPETAG_MASS>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
