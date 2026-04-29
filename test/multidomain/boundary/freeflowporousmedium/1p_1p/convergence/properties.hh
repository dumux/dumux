// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief The properties for a simple Darcy test (cell-centered finite volume method)
 */
#ifndef DUMUX_DARCYSTOKES_PROPERTIES_HH
#define DUMUX_DARCYSTOKES_PROPERTIES_HH

#ifndef TYPETAG_FFMOMENTUM
#define TYPETAG_FFMOMENTUM FreeFlowOnePMomentum
#endif

#ifndef TYPETAG_FFMASS
#define TYPETAG_FFMASS FreeFlowOnePMass
#endif

#ifndef TYPETAG_PMMASS
#define TYPETAG_PMMASS DarcyOneP
#endif

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/variables.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/hybrid/gridvariablescache.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/gridvariables.hh>

#include <dumux/multidomain/boundary/freeflowporousmedium/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "spatialparams.hh"
#include "problem_darcy.hh"
#include "problem_freeflow.hh"
#include "problem_freeflow_newinterface.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DarcyOneP { using InheritsFrom = std::tuple<OneP, CCTpfaModel>; };
struct DarcyOnePBox { using InheritsFrom = std::tuple<OneP, BoxModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_PMMASS>
{ using type = DarcySubProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::TYPETAG_PMMASS>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::TYPETAG_PMMASS>
{
#ifndef GRIDTYPE
    using type = Dune::YaspGrid<2>;
#else
    using type = GRIDTYPE;
#endif
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::TYPETAG_PMMASS>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ConvergenceTestSpatialParams<GridGeometry, Scalar>;
};

// Create new type tags
namespace TTag {
struct FreeFlowOneP {};
struct FreeFlowOnePMass { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMassOneP, CCTpfaModel>; };
struct FreeFlowOnePMomentum { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct FreeFlowOnePMassBox { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMassOneP, BoxModel>; };
struct FreeFlowOnePMomentumPQ1Bubble { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMomentumCVFE, PQ1BubbleModel>; };
struct FreeFlowOnePMomentumPQ1BubbleHybrid { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMomentumCVFE, PQ1BubbleHybridModel>; };
struct FreeFlowOnePMomentumPQ2Hybrid { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMomentumCVFE, PQ2HybridModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::FreeFlowOneP>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::FreeFlowOneP>
{
#ifndef GRIDTYPE
    using Coords = Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2>;
    using type = Dune::YaspGrid<2, Coords>;
#else
    using type = GRIDTYPE;
#endif
};

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMomentum>
{ using type = FreeFlowSubProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMass>
{ using type = FreeFlowSubProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMomentumPQ1Bubble>
{ using type = FreeFlowSubProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMomentumPQ1BubbleHybrid>
{ using type = FreeFlowSubProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMomentumPQ2Hybrid>
{ using type = FreeFlowSubProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMassBox>
{ using type = FreeFlowSubProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMassProblem<TypeTag>>; };

// the variables
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::FreeFlowOnePMomentumPQ1BubbleHybrid>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesMomentumCVFEVolumeVariablesTraits<PV, FSY, FST, MT>;
public:
    using type = NavierStokesMomentumCVFEVariables<Traits>;
};

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::FreeFlowOnePMomentumPQ2Hybrid>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;

    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    static_assert(!FSY::isMiscible(), "The Navier-Stokes model only works with immiscible fluid systems.");

    using Traits = NavierStokesMomentumCVFEVolumeVariablesTraits<PV, FSY, FST, MT>;
public:
    using type = NavierStokesMomentumCVFEVariables<Traits>;
};

//! The grid variables
template<class TypeTag>
struct GridVariables<TypeTag, TTag::FreeFlowOnePMomentumPQ1Bubble>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    // ToDo: Do not determine enableCache by EnableGridVolumeVariablesCache
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
struct GridVariables<TypeTag, TTag::FreeFlowOnePMomentumPQ1BubbleHybrid>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    // ToDo: Do not determine enableCache by EnableGridVolumeVariablesCache
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

template<class TypeTag>
struct GridVariables<TypeTag, TTag::FreeFlowOnePMomentumPQ2Hybrid>
{
private:
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    // ToDo: Do not determine enableCache by EnableGridVolumeVariablesCache
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridVolumeVariablesCache>();
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Variables = Dumux::Detail::CVFE::VariablesAdapter<GetPropType<TypeTag, Properties::VolumeVariables>>;
    using IPDataCache = Dumux::CVFE::LocalBasisInterpolationPointData<GG>;
    using Traits = Dumux::Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<Problem, Variables, IPDataCache>;
    using GVC = Dumux::Experimental::CVFE::HybridCVFEGridVariablesCache<Traits, enableCache>;
public:
    using type = Dumux::Experimental::GridVariables<GG, GVC>;
};

//! The grid variables
template<class TypeTag>
struct GridVariables<TypeTag, TTag::FreeFlowOnePMassBox>
{
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    // ToDo: Do not determine enableCache by EnableGridVolumeVariablesCache
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
struct EnableGridGeometryCache<TypeTag, TTag::FreeFlowOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FreeFlowOneP> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FreeFlowOneP> { static constexpr bool value = true; };

#ifdef QUADRATURE_RULE_FFMASS
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::FreeFlowOnePMassBox>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using QuadTraits = BoxQuadratureTraits<GridView, QUADRATURE_RULE_FFMASS, QUADRATURE_RULE_FFMASS>;

public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, BoxDefaultGridGeometryTraits<GridView, DefaultMapperTraits<GridView>, QuadTraits>>;
};
#endif

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TYPETAG_PMMASS>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_FFMOMENTUM, TTag::TYPETAG_FFMASS, TTag::TYPETAG_PMMASS>;
    using type = FreeFlowPorousMediumCouplingManager<Traits>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::FreeFlowOneP>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_FFMOMENTUM, TTag::TYPETAG_FFMASS, TTag::TYPETAG_PMMASS>;
    using type = FreeFlowPorousMediumCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
