// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Partition of a network following a distributed bulk grid, checked against a brute-force
 *        reference: an element belongs to the rank whose interior bulk elements contain its center,
 *        otherwise to the rank with the closest bulk domain boundary.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/foamgrid/foamgrid.hh>
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/common/initialize.hh>
#include <dumux/geometry/distance.hh>
#include <dumux/multidomain/embedded/lowdimpartition.hh>

namespace {

using Network = Dune::FoamGrid<1, 3>;
using Point = Dune::FieldVector<double, 3>;

//! random segments in [-0.3, 1.3]^3 and a line through the unit cube whose element centers lie on cell faces
std::unique_ptr<Network> makeNetwork()
{
    Dune::GridFactory<Network> factory;
    std::mt19937 generator(42);
    std::uniform_real_distribution<double> coordinate(-0.3, 1.3);
    unsigned int v = 0;
    for (int s = 0; s < 400; ++s)
    {
        factory.insertVertex(Point{coordinate(generator), coordinate(generator), coordinate(generator)});
        factory.insertVertex(Point{coordinate(generator), coordinate(generator), coordinate(generator)});
        factory.insertElement(Dune::GeometryTypes::line, {v, v + 1});
        v += 2;
    }
    const int n = 16;
    for (int i = 0; i <= n; ++i)
        factory.insertVertex(Point{0.5, 0.25, -0.25 + 1.5*i/n});
    for (int i = 0; i < n; ++i)
        factory.insertElement(Dune::GeometryTypes::line, {v + i, v + i + 1});
    return factory.createGrid();
}

//! squared distance of a point to an axis-aligned box given by its lower and upper corner
double squaredDistanceToBox(const Point& p, const Point& lo, const Point& hi)
{
    double d = 0.0;
    for (int k = 0; k < 3; ++k)
    {
        const double outside = std::max({lo[k] - p[k], 0.0, p[k] - hi[k]});
        d += outside*outside;
    }
    return d;
}

//! brute-force reference of this rank's claim per network element: 0 if an interior bulk element
//! contains the center, else the squared distance to the closest domain boundary facet of an
//! interior bulk element
template<class BulkGridView, class Contains, class FacetDistance>
std::vector<double> referenceClaims(const BulkGridView& bulkGridView, const Network& network,
                                    const Contains& contains, const FacetDistance& facetDistance)
{
    const auto networkView = network.leafGridView();
    std::vector<double> claim(networkView.size(0), std::numeric_limits<double>::max());
    for (const auto& e : elements(networkView))
    {
        const auto c = e.geometry().center();
        auto& d = claim[networkView.indexSet().index(e)];
        for (const auto& element : elements(bulkGridView, Dune::Partitions::interior))
        {
            if (contains(element, c))
            {
                d = 0.0;
                break;
            }
            for (const auto& intersection : intersections(bulkGridView, element))
                if (intersection.boundary())
                    d = std::min(d, facetDistance(intersection.geometry(), c));
        }
    }
    return claim;
}

template<class BulkGridView, class Contains, class FacetDistance>
int check(const std::string& name, const BulkGridView& bulkGridView, const Network& network,
          const Contains& contains, const FacetDistance& facetDistance)
{
    const auto& comm = bulkGridView.comm();
    const auto owner = Dumux::lowDimPartition(bulkGridView, network.leafGridView());
    const auto claim = referenceClaims(bulkGridView, network, contains, facetDistance);
    auto minClaim = claim;
    comm.min(minClaim.data(), minClaim.size());

    int failures = 0;
    std::uint64_t checksum = 0;
    int inside = 0, outside = 0;
    for (std::size_t i = 0; i < owner.size(); ++i)
    {
        checksum = checksum*1000003u + std::uint64_t(owner[i]) + 1u;
        if (owner[i] < 0 || owner[i] >= comm.size())
            ++failures;
        if (owner[i] == comm.rank() && claim[i] > minClaim[i] + 1e-12*(1.0 + minClaim[i]))
        {
            std::cerr << name << ": rank " << comm.rank() << " owns element " << i << " with squared distance "
                      << claim[i] << " but the closest is " << minClaim[i] << "\n";
            ++failures;
        }
        (minClaim[i] == 0.0 ? inside : outside) += 1;
    }
    if (comm.min(checksum) != comm.max(checksum))
        ++failures;

    std::vector<int> owned(comm.size(), 0);
    for (const auto r : owner)
        ++owned[r];

    failures = comm.sum(failures);
    if (comm.rank() == 0)
    {
        std::cout << name << ": " << inside << " elements inside, " << outside << " outside; per rank";
        for (const auto n : owned)
            std::cout << " " << n;
        std::cout << "; " << failures << " failure(s)\n";
    }
    return failures + (inside == 0) + (outside == 0);
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    const auto& helper = Dune::MPIHelper::instance();
    const auto network = makeNetwork();
    int failures = 0;

    {
        using Bulk = Dune::YaspGrid<3>;
        Bulk bulk(Point(1.0), std::array<int, 3>{{8, 8, 8}}, std::bitset<3>(), 1, helper.getCommunicator());
        const auto contains = [](const auto& element, const Point& p){
            const auto g = element.geometry();
            return squaredDistanceToBox(p, g.corner(0), g.corner(7)) == 0.0;
        };
        const auto facetDistance = [](const auto& facet, const Point& p){
            return squaredDistanceToBox(p, facet.corner(0), facet.corner(3));
        };
        failures += check("YaspGrid", bulk.leafGridView(), *network, contains, facetDistance);
    }

#if HAVE_DUNE_ALUGRID
    {
        using Bulk = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
        auto bulk = Dune::StructuredGridFactory<Bulk>::createSimplexGrid(Point(0.0), Point(1.0), std::array<unsigned int, 3>{{6, 6, 6}});
        bulk->loadBalance();
        const auto contains = [](const auto& element, const Point& p){
            const auto g = element.geometry();
            const auto local = g.local(p);
            const auto ref = Dune::referenceElement(g);
            return ref.checkInside(local);
        };
        const auto facetDistance = [](const auto& facet, const Point& p){
            return Dumux::squaredDistancePointTriangle(p, facet.corner(0), facet.corner(1), facet.corner(2));
        };
        failures += check("ALUGrid", bulk->leafGridView(), *network, contains, facetDistance);
    }
#endif

    return failures == 0 ? 0 : 1;
}
