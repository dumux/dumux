// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Test that the default grid variables select the cache depending on the element discretization
 *        and the variables depending on the type tag
 */
#include <config.h>

#include <type_traits>

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/properties.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/fcdiamond.hh>
#include <dumux/discretization/pq1.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/pq1nonconforming.hh>
#include <dumux/discretization/pq2.hh>
#include <dumux/discretization/pq3.hh>
#include <dumux/discretization/defaultgridvariables.hh>

namespace Dumux {

struct TestProblem {};
struct TestVolumeVariables {};
struct TestVariables {};

} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
struct DefaultGridVariablesTest {};
struct TestBox { using InheritsFrom = std::tuple<DefaultGridVariablesTest, BoxModel>; };
struct TestPQ1FE { using InheritsFrom = std::tuple<DefaultGridVariablesTest, PQ1FEModel>; };
struct TestPQ1Bubble { using InheritsFrom = std::tuple<DefaultGridVariablesTest, PQ1BubbleModel>; };
struct TestPQ1BubbleHybrid { using InheritsFrom = std::tuple<DefaultGridVariablesTest, PQ1BubbleHybridModel>; };
struct TestPQ1BubbleFE { using InheritsFrom = std::tuple<DefaultGridVariablesTest, PQ1BubbleFEModel>; };
struct TestDiamond { using InheritsFrom = std::tuple<DefaultGridVariablesTest, FaceCenteredDiamondModel>; };
struct TestPQ1NonconformingFE { using InheritsFrom = std::tuple<DefaultGridVariablesTest, PQ1NonconformingFEModel>; };
struct TestPQ2Hybrid { using InheritsFrom = std::tuple<DefaultGridVariablesTest, PQ2HybridModel>; };
struct TestPQ2FE { using InheritsFrom = std::tuple<DefaultGridVariablesTest, PQ2FEModel>; };
struct TestPQ3Hybrid { using InheritsFrom = std::tuple<DefaultGridVariablesTest, PQ3HybridModel>; };
struct TestPQ3FE { using InheritsFrom = std::tuple<DefaultGridVariablesTest, PQ3FEModel>; };

struct TestModelVariables {};
struct TestBoxUseVariables { using InheritsFrom = std::tuple<TestBox>; };
struct TestBoxWithVariables { using InheritsFrom = std::tuple<TestModelVariables, TestBoxUseVariables>; };
struct TestDiamondUseVariables { using InheritsFrom = std::tuple<TestDiamond>; };
struct TestPQ1BubbleUseVariables { using InheritsFrom = std::tuple<TestPQ1Bubble>; };
struct TestPQ2HybridWithVariables { using InheritsFrom = std::tuple<TestModelVariables, TestPQ2Hybrid>; };
struct TestPQ2FEWithVariables { using InheritsFrom = std::tuple<TestModelVariables, TestPQ2FE>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::DefaultGridVariablesTest> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::DefaultGridVariablesTest> { using type = double; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DefaultGridVariablesTest> { static constexpr bool value = false; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DefaultGridVariablesTest> { using type = TestProblem; };

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::DefaultGridVariablesTest> { using type = TestVolumeVariables; };

template<class TypeTag>
struct Variables<TypeTag, TTag::TestModelVariables> { using type = TestVariables; };

template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::TestBoxUseVariables> { static constexpr bool value = true; };

template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::TestDiamondUseVariables> { static constexpr bool value = true; };

template<class TypeTag>
struct Transition::UseVariables<TypeTag, TTag::TestPQ1BubbleUseVariables> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

namespace Dumux {

enum class CacheKind { controlVolume, hybrid, finiteElement };

template<class GG, CacheKind kind, class Variables>
using ExpectedCache = std::conditional_t<kind == CacheKind::controlVolume,
    Experimental::CVFE::CVFEGridVariablesCache<
        Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<TestProblem, Variables, CVFE::LocalBasisInterpolationPointData<GG>>, false>,
    std::conditional_t<kind == CacheKind::hybrid,
        Experimental::CVFE::HybridCVFEGridVariablesCache<
            Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<TestProblem, Variables, CVFE::LocalBasisInterpolationPointData<GG>>, false>,
        Experimental::FEGridVariablesCache<
            Experimental::FEDefaultGridVariablesCacheTraits<TestProblem, Variables, CVFE::LocalBasisInterpolationPointData<GG>>, false>
    >
>;

//! The cache selected for the given variables
template<class TypeTag, CacheKind kind>
constexpr bool picks()
{
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using Cache = Experimental::DefaultGridVariablesCache<GG, TestProblem, TestVariables, false>;
    return std::is_same_v<Cache, ExpectedCache<GG, kind, TestVariables>>;
}

//! The grid variables property of the type tag assembles with the given variables per local dof
template<class TypeTag, CacheKind kind, class Variables>
constexpr bool assemblesWith()
{
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using Expected = Experimental::GridVariables<GG, ExpectedCache<GG, kind, Variables>>;
    return Detail::usesVariables<TypeTag>()
        && std::is_same_v<GetPropType<TypeTag, Properties::GridVariables>, Expected>;
}

//! The grid variables property of the type tag assembles with the volume variables
template<class TypeTag>
constexpr bool assemblesWithVolumeVariables()
{
    using Expected = FVGridVariables<GetPropType<TypeTag, Properties::GridGeometry>,
                                     GetPropType<TypeTag, Properties::GridVolumeVariables>,
                                     GetPropType<TypeTag, Properties::GridFluxVariablesCache>>;
    return !Detail::usesVariables<TypeTag>()
        && std::is_same_v<GetPropType<TypeTag, Properties::GridVariables>, Expected>;
}

using namespace Properties::TTag;
using AdaptedVolumeVariables = Detail::CVFE::VariablesAdapter<TestVolumeVariables>;

static_assert(picks<TestBox, CacheKind::controlVolume>());
static_assert(picks<TestPQ1Bubble, CacheKind::controlVolume>());
static_assert(picks<TestDiamond, CacheKind::controlVolume>());
static_assert(picks<TestPQ1BubbleHybrid, CacheKind::hybrid>());
static_assert(picks<TestPQ2Hybrid, CacheKind::hybrid>());
static_assert(picks<TestPQ3Hybrid, CacheKind::hybrid>());
static_assert(picks<TestPQ1FE, CacheKind::finiteElement>());
static_assert(picks<TestPQ1BubbleFE, CacheKind::finiteElement>());
static_assert(picks<TestPQ1NonconformingFE, CacheKind::finiteElement>());
static_assert(picks<TestPQ2FE, CacheKind::finiteElement>());
static_assert(picks<TestPQ3FE, CacheKind::finiteElement>());

// schemes with a sub-control volume for every local dof keep the volume variables unless the flag is set
static_assert(assemblesWithVolumeVariables<TestBox>());
static_assert(assemblesWithVolumeVariables<TestPQ1Bubble>());
static_assert(assemblesWithVolumeVariables<TestDiamond>());
static_assert(assemblesWith<TestBoxUseVariables, CacheKind::controlVolume, AdaptedVolumeVariables>());
static_assert(assemblesWith<TestPQ1BubbleUseVariables, CacheKind::controlVolume, AdaptedVolumeVariables>());
static_assert(assemblesWith<TestDiamondUseVariables, CacheKind::controlVolume, AdaptedVolumeVariables>());
static_assert(assemblesWith<TestBoxWithVariables, CacheKind::controlVolume, TestVariables>());

// the other schemes always assemble with the variables per local dof, the model's if it defines them
static_assert(assemblesWith<TestPQ2Hybrid, CacheKind::hybrid, AdaptedVolumeVariables>());
static_assert(assemblesWith<TestPQ2HybridWithVariables, CacheKind::hybrid, TestVariables>());
static_assert(assemblesWith<TestPQ1FE, CacheKind::finiteElement, AdaptedVolumeVariables>());
static_assert(assemblesWith<TestPQ2FEWithVariables, CacheKind::finiteElement, TestVariables>());

} // end namespace Dumux

int main()
{
    return 0;
}
