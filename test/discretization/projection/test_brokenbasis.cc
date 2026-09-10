// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test the broken Lagrange basis over a topology-less entity set.
 */
#include <config.h>

#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/concepts/functionspacebasis_.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/discretization/projection/brokenbasis.hh>

namespace {

using Segment = Dune::MultiLinearGeometry<double, 1, 2>;
using SegmentSet = Dumux::GeometriesEntitySet<Segment>;

//! Two segments tiling the unit interval along the x-axis
auto makeSegmentSet()
{
    using P = Dune::FieldVector<double, 2>;
    std::vector<Segment> segments;
    segments.emplace_back(Dune::GeometryTypes::line, std::vector<P>{{0.0, 0.0}, {0.5, 0.0}});
    segments.emplace_back(Dune::GeometryTypes::line, std::vector<P>{{0.5, 0.0}, {1.0, 0.0}});
    return std::make_shared<SegmentSet>(std::move(segments));
}

void checkPartitionOfUnity(const auto& basis)
{
    auto localView = basis.localView();
    for (const auto& entity : basis.entitySet())
    {
        localView.bind(entity);
        const auto& localBasis = localView.tree().finiteElement().localBasis();
        std::vector<Dune::FieldVector<double, 1>> values;
        for (const double x : {0.0, 0.17, 0.37, 0.5, 0.83, 1.0})
        {
            localBasis.evaluateFunction({x}, values);
            const auto sum = std::accumulate(
                values.begin(), values.end(), 0.0,
                [] (double a, const auto& v) { return a + v[0]; }
            );
            if (Dune::FloatCmp::ne(sum, 1.0, 1e-14))
                DUNE_THROW(Dune::MathError, "Partition of unity violated: " << sum << " at x = " << x);
        }
    }
}

void checkIndexMapIsBijective(const auto& basis)
{
    std::vector<std::size_t> seen;
    auto localView = basis.localView();
    for (const auto& entity : basis.entitySet())
    {
        localView.bind(entity);
        const auto numLocalDofs = localView.tree().finiteElement().localBasis().size();
        for (std::size_t i = 0; i < numLocalDofs; ++i)
            seen.push_back(localView.index(i));
    }

    if (seen.size() != basis.size())
        DUNE_THROW(Dune::InvalidStateException,
                   "Local dofs sum to " << seen.size() << " but basis reports " << basis.size());

    std::sort(seen.begin(), seen.end());
    std::vector<std::size_t> expected(basis.size());
    std::iota(expected.begin(), expected.end(), 0);
    if (seen != expected)
        DUNE_THROW(Dune::InvalidStateException, "Local-to-global map is not a bijection onto [0, size())");
}

//! On a broken basis no global index may be produced by more than one entity
void checkNoSharedDofs(const auto& basis)
{
    auto localView = basis.localView();
    std::vector<int> owners(basis.size(), -1);
    for (const auto& entity : basis.entitySet())
    {
        localView.bind(entity);
        const auto eIdx = static_cast<int>(basis.entitySet().index(entity));
        const auto numLocalDofs = localView.tree().finiteElement().localBasis().size();
        for (std::size_t i = 0; i < numLocalDofs; ++i)
        {
            const auto global = localView.index(i);
            if (owners[global] != -1 && owners[global] != eIdx)
                DUNE_THROW(Dune::InvalidStateException, "Dof " << global << " is shared between entities");
            owners[global] = eIdx;
        }
    }
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);

    using Basis1 = BrokenLagrangeBasis<SegmentSet, 1>;
    using Basis2 = BrokenLagrangeBasis<SegmentSet, 2>;

    static_assert(Concept::ProjectionBasis<Basis1>);
    static_assert(Concept::ProjectionBasis<Basis2>);
    static_assert(Concept::EntityRangeProvider<Basis2>);
    static_assert(Concept::ScalarLocalFiniteElement<typename Basis2::FiniteElement>);
    static_assert(Basis2::dimension == 1);

    // a type missing localView() must not satisfy the concept
    struct NotABasis { std::size_t size() const { return 0; } };
    static_assert(!Concept::ProjectionBasis<NotABasis>);
    static_assert(!Concept::EntityRangeProvider<NotABasis>);

    const auto segments = makeSegmentSet();

    const Basis1 basis1{segments};
    if (basis1.size() != 4)
        DUNE_THROW(Dune::InvalidStateException, "Expected 4 dofs for P1 on two segments, got " << basis1.size());

    const Basis2 basis2{segments};
    if (basis2.size() != 6)
        DUNE_THROW(Dune::InvalidStateException, "Expected 6 dofs for P2 on two segments, got " << basis2.size());

    checkPartitionOfUnity(basis1);
    checkPartitionOfUnity(basis2);
    checkIndexMapIsBijective(basis1);
    checkIndexMapIsBijective(basis2);
    checkNoSharedDofs(basis1);
    checkNoSharedDofs(basis2);

    // mixed geometry types in one set: a triangle and a quadrilateral in 3d
    {
        using Facet = Dune::MultiLinearGeometry<double, 2, 3>;
        using P = Dune::FieldVector<double, 3>;
        std::vector<Facet> facets;
        facets.emplace_back(Dune::GeometryTypes::triangle,
                            std::vector<P>{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}});
        facets.emplace_back(Dune::GeometryTypes::quadrilateral,
                            std::vector<P>{{1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {2.0, 1.0, 0.0}});
        using FacetSet = GeometriesEntitySet<Facet>;
        const auto facetSet = std::make_shared<FacetSet>(std::move(facets));

        const BrokenLagrangeBasis<FacetSet, 2> mixed{facetSet};
        // P2 on a triangle has 6 dofs, Q2 on a quadrilateral has 9
        if (mixed.size() != 15)
            DUNE_THROW(Dune::InvalidStateException, "Expected 15 dofs for mixed P2/Q2, got " << mixed.size());

        checkPartitionOfUnity(mixed);
        checkIndexMapIsBijective(mixed);
        checkNoSharedDofs(mixed);
    }

    std::cout << "All broken basis checks passed" << std::endl;
    return 0;
}
