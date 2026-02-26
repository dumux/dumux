// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid (Navier-)Stokes model with analytical solution (Donea 2003, \cite Donea2003).
 */
#ifndef DUMUX_DONEA_TEST_PROPERTIES_MOMENTUM_HH
#define DUMUX_DONEA_TEST_PROPERTIES_MOMENTUM_HH

#ifndef ENABLECACHING
#define ENABLECACHING true
#endif

#ifndef ENABLEFLUXVARSCACHING
#define ENABLEFLUXVARSCACHING ENABLECACHING
#endif

#ifndef GRIDTYPE
#define GRIDTYPE Dune::YaspGrid<2>
#endif

#ifndef NAVIER_STOKES_MODEL
#define NAVIER_STOKES_MODEL NavierStokesMomentum
#endif

#ifndef DISCRETIZATION_MODEL
#define DISCRETIZATION_MODEL FaceCenteredStaggeredModel
#endif

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM DoneaTestMomentum
#endif

#ifndef NEW_PROBLEM_INTERFACE
#define NEW_PROBLEM_INTERFACE 0
#endif

#ifndef NEW_VARIABLES_INTERFACE
#define NEW_VARIABLES_INTERFACE 0
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/flux/fluxvariablescaching.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/variables.hh>

#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/discretization/pq1bubble/fvelementgeometry.hh>
#include <dumux/discretization/pq2.hh>

#include "problem.hh"
#include "problem_newinterface.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DoneaTest { };
struct DoneaTestMomentum { using InheritsFrom = std::tuple<DoneaTest, NAVIER_STOKES_MODEL, DISCRETIZATION_MODEL>; };
struct DoneaTestMomentumPQ1Bubble { using InheritsFrom = std::tuple<DoneaTest, NavierStokesMomentumCVFE, PQ1BubbleModel>; };
struct DoneaTestMomentumPQ1BubbleHybrid { using InheritsFrom = std::tuple<DoneaTest, NavierStokesMomentumCVFE, PQ1BubbleHybridModel>; };
struct DoneaTestMomentumPQ2Hybrid { using InheritsFrom = std::tuple<DoneaTest, NavierStokesMomentumCVFE, PQ2HybridModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{
#if NEW_PROBLEM_INTERFACE
    using type = Dumux::DoneaTestProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>;
#else
    using type = Dumux::DoneaTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>;
#endif
};

// the fluid system
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

#ifdef QUADRATURE_RULE
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
    using QuadTraits = PQ1BubbleQuadratureTraits<GridView, QUADRATURE_RULE, QUADRATURE_RULE>;
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
                                           QUADRATURE_RULE,
                                           QUADRATURE_RULE,
                                           QUADRATURE_RULE>;

public:
    using type = PQ2FVGridGeometry<Scalar, GridView, enableCache, PQ2DefaultGridGeometryTraits<GridView, PQ2MapperTraits<GridView>, QuadTraits>>;
};
#endif

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DoneaTest> { using type = GRIDTYPE; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TYPETAG_MOMENTUM> { static constexpr bool value = ENABLECACHING; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TYPETAG_MOMENTUM> { static constexpr bool value = ENABLEFLUXVARSCACHING; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TYPETAG_MOMENTUM> { static constexpr bool value = ENABLECACHING; };

} // end namespace Dumux::Properties

#endif
