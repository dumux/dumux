// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief The properties for the coupled free-flow/porous-medium flow convergence test
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
#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq1.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/pq3.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include <dumux/multidomain/boundary/freeflowporousmedium/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "spatialparams.hh"
#include "problem_darcy.hh"
#include "problem_darcy_newinterface.hh"
#include "problem_freeflow.hh"
#include "problem_freeflow_newinterface.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DarcyOnePCommon { using InheritsFrom = std::tuple<OneP>; };
struct DarcyOneP { using InheritsFrom = std::tuple<DarcyOnePCommon, CCTpfaModel>; };

//! The porous medium domain assembled with the new interfaces
struct DarcyOnePBoxNewInterface { using InheritsFrom = std::tuple<DarcyOnePCommon, BoxModel>; };
struct DarcyOnePPQ2Hybrid { using InheritsFrom = std::tuple<DarcyOnePCommon, PQ2HybridModel>; };
struct DarcyOnePPQ3Hybrid { using InheritsFrom = std::tuple<DarcyOnePCommon, PQ3HybridModel>; };
struct DarcyOnePPQ1FE { using InheritsFrom = std::tuple<DarcyOnePCommon, PQ1FEModel>; };
struct DarcyOnePPQ3FE { using InheritsFrom = std::tuple<DarcyOnePCommon, PQ3FEModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOnePCommon>
{ using type = DarcySubProblem<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyOnePCommon>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> > ;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DarcyOnePCommon>
{
#ifndef GRIDTYPE
    using type = Dune::YaspGrid<2>;
#else
    using type = GRIDTYPE;
#endif
};

template<class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyOnePCommon>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = ConvergenceTestSpatialParams<GridGeometry, Scalar>;
};

//! The problem stated in terms of the new interfaces
template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOnePBoxNewInterface>
{ using type = DarcySubProblemNewInterface<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOnePPQ2Hybrid>
{ using type = DarcySubProblemNewInterface<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOnePPQ3Hybrid>
{ using type = DarcySubProblemNewInterface<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOnePPQ1FE>
{ using type = DarcySubProblemNewInterface<TypeTag>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DarcyOnePPQ3FE>
{ using type = DarcySubProblemNewInterface<TypeTag>; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DarcyOnePPQ1FE> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DarcyOnePPQ3FE> { static constexpr bool value = true; };

//! Assemble with the variables per local dof
template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::DarcyOnePBoxNewInterface> { static constexpr bool value = true; };

// Create new type tags
namespace TTag {
struct FreeFlowOneP {};
struct FreeFlowOnePMass { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMassOneP, CCTpfaModel>; };
struct FreeFlowOnePMomentum { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct FreeFlowOnePMassBox { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMassOneP, BoxModel>; };
struct FreeFlowOnePMomentumCVFE { using InheritsFrom = std::tuple<FreeFlowOneP>; };
struct FreeFlowOnePMomentumPQ1Bubble { using InheritsFrom = std::tuple<FreeFlowOnePMomentumCVFE, NavierStokesMomentumCVFE, PQ1BubbleModel>; };
struct FreeFlowOnePMomentumPQ1BubbleHybrid { using InheritsFrom = std::tuple<FreeFlowOnePMomentumCVFE, NavierStokesMomentumCVFE, PQ1BubbleHybridModel>; };
struct FreeFlowOnePMomentumPQ2Hybrid { using InheritsFrom = std::tuple<FreeFlowOnePMomentumCVFE, NavierStokesMomentumCVFE, PQ2HybridModel>; };
struct FreeFlowOnePMomentumPQ3Hybrid { using InheritsFrom = std::tuple<FreeFlowOnePMomentumCVFE, NavierStokesMomentumCVFE, PQ3HybridModel>; };
struct FreeFlowOnePMassPQ2Hybrid { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMassOneP, PQ2HybridModel>; };

//! The pure finite element pairs: P2/P1 (Taylor-Hood) and P3/P2
struct FreeFlowOnePMomentumPQ2FE { using InheritsFrom = std::tuple<FreeFlowOnePMomentumCVFE, NavierStokesMomentumFE, PQ2FEModel>; };
struct FreeFlowOnePMomentumPQ3FE { using InheritsFrom = std::tuple<FreeFlowOnePMomentumCVFE, NavierStokesMomentumFE, PQ3FEModel>; };
struct FreeFlowOnePMassPQ1FE { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMassOnePFE, PQ1FEModel>; };
struct FreeFlowOnePMassPQ2FE { using InheritsFrom = std::tuple<FreeFlowOneP, NavierStokesMassOnePFE, PQ2FEModel>; };
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

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMomentumPQ3Hybrid>
{ using type = FreeFlowSubProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMomentumPQ2FE>
{ using type = FreeFlowSubProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMomentumPQ3FE>
{ using type = FreeFlowSubProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMassPQ2Hybrid>
{ using type = FreeFlowSubProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMassPQ1FE>
{ using type = FreeFlowSubProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMassProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowOnePMassPQ2FE>
{ using type = FreeFlowSubProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMassProblem<TypeTag>>; };

//! Assemble with the variables per local dof
template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::FreeFlowOnePMomentumPQ1Bubble> { static constexpr bool value = true; };

//! Assemble with the variables per local dof
template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::FreeFlowOnePMassBox> { static constexpr bool value = true; };

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

/*!
 * \brief Quadrature rules of sufficient order for the PQ3 velocity and the PQ2 pressure spaces
 * \note The PQ2 free-flow mass domain integrates the fluxes of a PQ3 velocity.
 */
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::FreeFlowOnePMomentumPQ3Hybrid>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Rule = Dumux::QuadratureRules::DuneQuadrature<8>;
    using QuadTraits = PQ3QuadratureTraits<GridView, Rule, Rule, Rule, Rule, Rule>;
public:
    using type = PQ3FVGridGeometry<Scalar, GridView, enableCache,
                                   PQ3DefaultGridGeometryTraits<GridView, PQ3MapperTraits<GridView>, QuadTraits>>;
};

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::FreeFlowOnePMassPQ2Hybrid>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Rule = Dumux::QuadratureRules::DuneQuadrature<6>;
    using QuadTraits = PQ2QuadratureTraits<GridView, Rule, Rule, Rule, Rule, Rule>;
public:
    using type = PQ2FVGridGeometry<Scalar, GridView, enableCache,
                                   PQ2DefaultGridGeometryTraits<GridView, PQ2MapperTraits<GridView>, QuadTraits>>;
};

template<class TypeTag>
struct GridGeometry<TypeTag, TTag::DarcyOnePPQ2Hybrid>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Rule = Dumux::QuadratureRules::DuneQuadrature<6>;
    using QuadTraits = PQ2QuadratureTraits<GridView, Rule, Rule, Rule, Rule, Rule>;
public:
    using type = PQ2FVGridGeometry<Scalar, GridView, enableCache,
                                   PQ2DefaultGridGeometryTraits<GridView, PQ2MapperTraits<GridView>, QuadTraits>>;
};

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DarcyOnePCommon>
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
