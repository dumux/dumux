// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesNCTests
 * \brief The properties of the test for the Navier-Stokes models with analytical solution.
 */
#ifndef DUMUX_SINCOS_TEST_PROPERTIES_HH
#define DUMUX_SINCOS_TEST_PROPERTIES_HH

#ifndef TYPETAG_MOMENTUM
#define TYPETAG_MOMENTUM SincosTestMomentum
#endif

#ifndef TYPETAG_MASS
#define TYPETAG_MASS SincosTestMass
#endif

#ifndef NEW_PROBLEM_INTERFACE
#define NEW_PROBLEM_INTERFACE 0
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/momentum/fcstaggered/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/model.hh>
#include <dumux/freeflow/navierstokes/momentum/cvfe/variables.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
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
#include "problem_newinterface.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct SincosTest {};
struct SincosTestMomentum { using InheritsFrom = std::tuple<SincosTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct SincosTestMass { using InheritsFrom = std::tuple<SincosTest, NavierStokesMassOneP, CCTpfaModel>; };
struct SincosTestMomentumPQ1BubbleHybrid { using InheritsFrom = std::tuple<SincosTest, NavierStokesMomentumCVFE, PQ1BubbleHybridModel>; };
struct SincosTestMomentumPQ2Hybrid { using InheritsFrom = std::tuple<SincosTest, NavierStokesMomentumCVFE, PQ2HybridModel>; };
struct SincosTestMassBox { using InheritsFrom = std::tuple<SincosTest, NavierStokesMassOneP, BoxModel>; };
struct SincosTestMomentumOnly { using InheritsFrom = std::tuple<SincosTest>; };
struct SincosTestMomentumOnlyPQ1BubbleHybrid { using InheritsFrom = std::tuple<SincosTestMomentumOnly, NavierStokesMomentumCVFE, PQ1BubbleHybridModel>; };
struct SincosTestMomentumOnlyPQ2Hybrid { using InheritsFrom = std::tuple<SincosTestMomentumOnly, NavierStokesMomentumCVFE, PQ2HybridModel>; };
} // end namespace TTag

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::SincosTest>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::SincosTest> { using type = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<GetPropType<TypeTag, Properties::Scalar>, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MOMENTUM>
{
#if NEW_PROBLEM_INTERFACE
    using type = SincosTestProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>;
#else
    using type = SincosTestProblem<TypeTag, Dumux::NavierStokesMomentumProblem<TypeTag>>;
#endif
};

template<class TypeTag>
struct Problem<TypeTag, TTag::TYPETAG_MASS>
{
#if NEW_PROBLEM_INTERFACE
    using type = SincosTestProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMassProblem<TypeTag>>;
#else
    using type = SincosTestProblem<TypeTag, Dumux::NavierStokesMassProblem<TypeTag>>;
#endif
};

#if NEW_PROBLEM_INTERFACE
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

template<class TypeTag>
struct GridVariables<TypeTag, TTag::SincosTestMassBox>
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
#endif

// the momentum problem without coupling to a mass problem
template<class TypeTag>
struct Problem<TypeTag, TTag::SincosTestMomentumOnly>
{ using type = SincosTestProblemNewInterface<TypeTag, Dumux::CVFENavierStokesMomentumProblem<TypeTag>>; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::SincosTestMomentumOnly>
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

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::SincosTestMomentumOnly>
{
    struct EmptyCouplingManager {};
    using type = EmptyCouplingManager;
};

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::SincosTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::SincosTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::SincosTest> { static constexpr bool value = true; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::SincosTest>
{
    using Traits = MultiDomainTraits<TTag::TYPETAG_MOMENTUM, TTag::TYPETAG_MASS>;
    using type = FreeFlowCouplingManager<Traits>;
};

} // end namespace Dumux::Properties

#endif
