// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief The properties of the Taylor-Green vortex test.
 */
#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_TAYLOR_GREEN_VORTEX_PROPERTIES_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_TAYLOR_GREEN_VORTEX_PROPERTIES_HH

// the spatial dimension (note that a macro named DIM would clash with dune-uggrid)
#ifndef TAYLORGREEN_DIM
#define TAYLORGREEN_DIM 2
#endif

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM TaylorGreenTestMomentumPQ1BubbleHybrid
#endif

#ifndef TYPETAG_MASS
#define TYPETAG_MASS TaylorGreenTestMassBox
#endif

// Non-hybrid PQ1Bubble (one bubble per element) combined with Box pressure is only stable on simplices
#ifndef SIMPLEX_GRID
#define SIMPLEX_GRID 0
#endif

#include <dune/grid/yaspgrid.hh>
#if SIMPLEX_GRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/variables.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/pq3.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/gridvariables.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem_2d.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct TaylorGreenTest {};
struct TaylorGreenTestMomentumPQ1Bubble { using InheritsFrom = std::tuple<TaylorGreenTest, NavierStokesMomentumCVFE, PQ1BubbleModel>; };
struct TaylorGreenTestMomentumPQ1BubbleHybrid { using InheritsFrom = std::tuple<TaylorGreenTest, NavierStokesMomentumCVFE, PQ1BubbleHybridModel>; };
struct TaylorGreenTestMomentumPQ2Hybrid { using InheritsFrom = std::tuple<TaylorGreenTest, NavierStokesMomentumCVFE, PQ2HybridModel>; };
struct TaylorGreenTestMomentumPQ2HybridHOQuad { using InheritsFrom = std::tuple<TaylorGreenTestMomentumPQ2Hybrid>; };
struct TaylorGreenTestMomentumPQ3Hybrid { using InheritsFrom = std::tuple<TaylorGreenTest, NavierStokesMomentumCVFE, PQ3HybridModel>; };
struct TaylorGreenTestMassBox { using InheritsFrom = std::tuple<TaylorGreenTest, NavierStokesMassOneP, BoxModel>; };
struct TaylorGreenTestMassBoxHOQuad { using InheritsFrom = std::tuple<TaylorGreenTestMassBox>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TaylorGreenTest>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TaylorGreenTest>
#if SIMPLEX_GRID
{ using type = Dune::ALUGrid<TAYLORGREEN_DIM, TAYLORGREEN_DIM, Dune::simplex, Dune::nonconforming>; };
#else
{ using type = Dune::YaspGrid<TAYLORGREEN_DIM, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, TAYLORGREEN_DIM> >; };
#endif

// Set the problem property
template<class TypeTag, class BaseProblem>
using TaylorGreenProblem = TaylorGreenVortexProblem2D<TypeTag, BaseProblem>;

template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{ using type = TaylorGreenProblem<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TaylorGreenTestMassBox>
{ using type = TaylorGreenProblem<TypeTag, Dumux::CVFENavierStokesMassProblem<TypeTag>>; };

// the momentum variables
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TYPETAG_MOMENTUM>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    using Traits = NavierStokesMomentumCVFEVolumeVariablesTraits<PV, FSY, FST, MT>;
public:
    using type = NavierStokesMomentumCVFEVariables<Traits>;
};

// the grid variables of the non-hybrid schemes (hybrid schemes set appropriate defaults)
template<class TypeTag, class GG>
struct TaylorGreenCVFEGridVariables
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::CVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::TaylorGreenTestMassBox>
: public TaylorGreenCVFEGridVariables<TypeTag, GetPropType<TypeTag, Properties::GridGeometry>> {};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::TaylorGreenTestMomentumPQ1Bubble>
: public TaylorGreenCVFEGridVariables<TypeTag, GetPropType<TypeTag, Properties::GridGeometry>> {};

// Hybrid PQ2 with higher-order quadrature rules: by default, the source term is integrated over the
// sub-control volumes (of the vertex degrees of freedom) with the midpoint rule, which is only
// second-order consistent. Here, a Gauss rule is used instead (which does not change the
// convergence order of hybrid PQ2, see README.md).
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::TaylorGreenTestMomentumPQ2HybridHOQuad>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using QuadTraits = PQ2QuadratureTraits<GridView,
                                           QuadratureRules::DuneQuadrature<3>, // sub-control volumes
                                           QuadratureRules::DuneQuadrature<2>, // sub-control-volume faces
                                           QuadratureRules::DuneQuadrature<4>, // elements
                                           QuadratureRules::DuneQuadrature<4>, // intersections
                                           QuadratureRules::DuneQuadrature<4>>; // boundary faces
public:
    using type = PQ2FVGridGeometry<Scalar, GridView, enableCache,
                                   PQ2DefaultGridGeometryTraits<GridView, PQ2MapperTraits<GridView>, QuadTraits>>;
};

// Box scheme for the mass balance with a Gauss rule for the mass fluxes over the sub-control-volume
// faces, which integrates the (higher-order) velocity of the momentum discretization more accurately
// than the default midpoint rule
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::TaylorGreenTestMassBoxHOQuad>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using QuadTraits = BoxQuadratureTraits<GridView, QuadratureRules::MidpointQuadrature, QuadratureRules::DuneQuadrature<2>>;
public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache,
                                   BoxDefaultGridGeometryTraits<GridView, DefaultMapperTraits<GridView>, QuadTraits>>;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TaylorGreenTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TaylorGreenTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TaylorGreenTest> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TaylorGreenTest>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_MOMENTUM, TTag::TYPETAG_MASS>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
