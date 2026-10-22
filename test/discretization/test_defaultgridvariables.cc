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
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::DefaultGridVariablesTest> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::DefaultGridVariablesTest> { using type = double; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DefaultGridVariablesTest> { static constexpr bool value = false; };

} // end namespace Dumux::Properties

namespace Dumux {

struct TestProblem {};
struct TestVariables {};

enum class CacheKind { controlVolume, hybrid, finiteElement };

template<class TypeTag, CacheKind kind>
constexpr bool picks()
{
    using GG = GetPropType<TypeTag, Properties::GridGeometry>;
    using IPD = CVFE::LocalBasisInterpolationPointData<GG>;
    using Cache = Experimental::DefaultGridVariablesCache<GG, TestProblem, TestVariables, false>;

    using Expected = std::conditional_t<kind == CacheKind::controlVolume,
        Experimental::CVFE::CVFEGridVariablesCache<Experimental::CVFE::CVFEDefaultGridVariablesCacheTraits<TestProblem, TestVariables, IPD>, false>,
        std::conditional_t<kind == CacheKind::hybrid,
            Experimental::CVFE::HybridCVFEGridVariablesCache<Experimental::CVFE::HybridCVFEDefaultGridVariablesCacheTraits<TestProblem, TestVariables, IPD>, false>,
            Experimental::FEGridVariablesCache<Experimental::FEDefaultGridVariablesCacheTraits<TestProblem, TestVariables, IPD>, false>
        >
    >;

    return std::is_same_v<Cache, Expected>;
}

using namespace Properties::TTag;
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

} // end namespace Dumux

int main()
{
    return 0;
}
