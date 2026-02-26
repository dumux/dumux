// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties of the test for the staggered grid (Navier-)Stokes model with analytical solution (Donea 2003, \cite Donea2003).
 */
#ifndef DUMUX_DONEA_TEST_PROPERTIES_HH
#define DUMUX_DONEA_TEST_PROPERTIES_HH

#ifndef ENABLECACHING
#define ENABLECACHING 0
#endif

#ifndef ENABLEFLUXVARSCACHING
#define ENABLEFLUXVARSCACHING ENABLECACHING
#endif

#ifndef NAVIER_STOKES_MOMENTUM_MODEL
#define NAVIER_STOKES_MOMENTUM_MODEL NavierStokesMomentum
#endif

#ifndef MOMENTUM_DISCRETIZATION_MODEL
#define MOMENTUM_DISCRETIZATION_MODEL FaceCenteredStaggeredModel
#endif

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM DoneaTestMomentum
#endif

#ifndef TYPETAG_MASS
#define TYPETAG_MASS DoneaTestMass
#endif

#ifndef MASS_DISCRETIZATION_MODEL
#define MASS_DISCRETIZATION_MODEL CCTpfaModel
#endif

#ifndef ALUGRID_CELL_TYPE
#define ALUGRID_CELL_TYPE cube
#endif

#ifndef NEW_PROBLEM_INTERFACE
#define NEW_PROBLEM_INTERFACE 0
#endif

#ifndef NEW_VARIABLES_INTERFACE
#define NEW_VARIABLES_INTERFACE 0
#endif

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#else
#include <dune/grid/yaspgrid.hh>
#endif

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq2.hh>

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/variables.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"
#include "problem_newinterface.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DoneaTest {};
struct DoneaTestMomentum { using InheritsFrom = std::tuple<DoneaTest, NAVIER_STOKES_MOMENTUM_MODEL, MOMENTUM_DISCRETIZATION_MODEL>; };
struct DoneaTestMomentumPQ1Bubble { using InheritsFrom = std::tuple<DoneaTest, NavierStokesMomentumCVFE, PQ1BubbleModel>; };
struct DoneaTestMomentumPQ1BubbleHybrid { using InheritsFrom = std::tuple<DoneaTest, NavierStokesMomentumCVFE, PQ1BubbleHybridModel>; };
struct DoneaTestMomentumPQ2Hybrid { using InheritsFrom = std::tuple<DoneaTest, NavierStokesMomentumCVFE, PQ2HybridModel>; };
struct DoneaTestMass { using InheritsFrom = std::tuple<DoneaTest, NavierStokesMassOneP, MASS_DISCRETIZATION_MODEL>; };
struct DoneaTestMassBox { using InheritsFrom = std::tuple<DoneaTest, NavierStokesMassOneP, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{
#if NEW_PROBLEM_INTERFACE
    using type = Dumux::DoneaTestProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>;
#else
    using type = Dumux::DoneaTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>;
#endif
};

template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MASS>
{
#if NEW_PROBLEM_INTERFACE
    using type = DoneaTestProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMassProblem<TypeTag>>;
#else
    using type = DoneaTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>;
#endif

};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DoneaTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

#if NEW_VARIABLES_INTERFACE
// the variables
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::TYPETAG_MOMENTUM>
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
#endif

#ifdef QUADRATURE_RULE_MOMENTUM
// use custom quadrature rules for hybrid pq1bubble scheme
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::DoneaTestMomentumPQ1BubbleHybrid>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    // Higher order quadrature rule only works for hybrid pq1bubble scheme
    // since we can't generate geometries for overlapping scvs
    using QuadTraits = PQ1BubbleQuadratureTraits<GridView, QUADRATURE_RULE_MOMENTUM, QUADRATURE_RULE_MOMENTUM>;
    using QuadratureGridGeometryTraits = HybridPQ1BubbleCVFEGridGeometryTraits<PQ1BubbleDefaultGridGeometryTraits<GridView, PQ1BubbleMapperTraits<GridView>, QuadTraits>>;

public:
    using type = PQ1BubbleFVGridGeometry<Scalar, GridView, enableCache, QuadratureGridGeometryTraits>;
};

// use custom quadrature rules for hybrid pq2 scheme
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::DoneaTestMomentumPQ2Hybrid>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using QuadTraits = PQ2QuadratureTraits<GridView,
                                           Dumux::QuadratureRules::MidpointQuadrature,
                                           QUADRATURE_RULE_MOMENTUM,
                                           QUADRATURE_RULE_MOMENTUM,
                                           QUADRATURE_RULE_MOMENTUM>;

public:
    using type = PQ2FVGridGeometry<Scalar, GridView, enableCache, PQ2DefaultGridGeometryTraits<GridView, PQ2MapperTraits<GridView>, QuadTraits>>;
};
#endif

#ifdef QUADRATURE_RULE_MASS
// use custom quadrature rules
template<class TypeTag>
struct GridGeometry<TypeTag, TTag::TYPETAG_MASS>
{
private:
    static constexpr bool enableCache = getPropValue<TypeTag, Properties::EnableGridGeometryCache>();
    using GridView = typename GetPropType<TypeTag, Properties::Grid>::LeafGridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using QuadTraits = BoxQuadratureTraits<GridView, Dumux::QuadratureRules::MidpointQuadrature, QUADRATURE_RULE_MASS>;

public:
    using type = BoxFVGridGeometry<Scalar, GridView, enableCache, BoxDefaultGridGeometryTraits<GridView, DefaultMapperTraits<GridView>, QuadTraits>>;
};
#endif


template<class TypeTag>
struct Grid<TypeTag, TTag::DoneaTest>
#if HAVE_DUNE_ALUGRID
{ using type = Dune::ALUGrid<2, 2, Dune::ALUGRID_CELL_TYPE, Dune::nonconforming>; };
#else
{ using type = Dune::YaspGrid<2>; };
#endif

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DoneaTest> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DoneaTest> { static constexpr bool value = ENABLEFLUXVARSCACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DoneaTest> { static constexpr bool value = ENABLECACHING; };

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DoneaTest>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_MOMENTUM, TTag::TYPETAG_MASS>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
