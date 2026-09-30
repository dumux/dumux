// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief The properties of the Taylor-Green vortex test.
 */
#ifndef DUMUX_TAYLOR_GREEN_VORTEX_TEST_PROPERTIES_HH
#define DUMUX_TAYLOR_GREEN_VORTEX_TEST_PROPERTIES_HH

// the spatial dimension (note that a macro named DIM would clash with dune-uggrid)
#ifndef TAYLORGREEN_DIM
#define TAYLORGREEN_DIM 2
#endif

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM TaylorGreenTestMomentumPQ1BubbleHybrid
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
#include <dumux/discretization/cvfe/gridvariablescache_.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/gridvariables.hh>

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct TaylorGreenTest {};
struct TaylorGreenTestMomentumPQ1Bubble { using InheritsFrom = std::tuple<TaylorGreenTest, NavierStokesMomentumCVFE, PQ1BubbleModel>; };
struct TaylorGreenTestMomentumPQ1BubbleHybrid { using InheritsFrom = std::tuple<TaylorGreenTest, NavierStokesMomentumCVFE, PQ1BubbleHybridModel>; };
struct TaylorGreenTestMomentumPQ2Hybrid { using InheritsFrom = std::tuple<TaylorGreenTest, NavierStokesMomentumCVFE, PQ2HybridModel>; };
struct TaylorGreenTestMassBox { using InheritsFrom = std::tuple<TaylorGreenTest, NavierStokesMassOneP, BoxModel>; };
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
template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{ using type = TaylorGreenTestProblem<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::TaylorGreenTestMassBox>
{ using type = TaylorGreenTestProblem<TypeTag, Dumux::CVFENavierStokesMassProblem<TypeTag>>; };

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

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::TaylorGreenTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::TaylorGreenTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::TaylorGreenTest> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::TaylorGreenTest>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_MOMENTUM, TTag::TaylorGreenTestMassBox>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
