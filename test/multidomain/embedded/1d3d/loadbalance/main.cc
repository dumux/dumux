// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Partitions of a bulk grid that balance the bulk elements together with the work of an embedded network
 *
 * The bulk grid is the cube [0,2]^3, the network lies in the octant [0,1]^3, so a partition of the bulk
 * elements alone gives the network to few processes. The bulk grid is created on rank 0, optionally refined,
 * and distributed with the load weights of the bulk elements, by ALUGrid's partitioner or by a Scotch
 * partition of the element graph. The sum of the weights over the processes has to equal the sum computed on
 * the undistributed grid, and the largest sum on a process may exceed the mean by at most a given factor.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/multidomain/embedded/bulkloadweights.hh>
#include <dumux/parallel/scotchpartitioner.hh>

namespace {

using Bulk = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>;
using Network = Dune::FoamGrid<1, 3>;
using Point = Dune::FieldVector<double, 3>;

//! the cube [0,2]^3 with n^3 hexahedra, all on rank 0
std::unique_ptr<Bulk> makeBulk(int n, int rank)
{
    Dune::GridFactory<Bulk> factory;
    if (rank == 0)
    {
        const auto vertex = [n](int i, int j, int k){ return static_cast<unsigned int>(i + (n+1)*(j + (n+1)*k)); };
        for (int k = 0; k <= n; ++k)
            for (int j = 0; j <= n; ++j)
                for (int i = 0; i <= n; ++i)
                    factory.insertVertex(Point{2.0*i/n, 2.0*j/n, 2.0*k/n});
        for (int k = 0; k < n; ++k)
            for (int j = 0; j < n; ++j)
                for (int i = 0; i < n; ++i)
                    factory.insertElement(Dune::GeometryTypes::cube(3), {
                        vertex(i, j, k), vertex(i+1, j, k), vertex(i, j+1, k), vertex(i+1, j+1, k),
                        vertex(i, j, k+1), vertex(i+1, j, k+1), vertex(i, j+1, k+1), vertex(i+1, j+1, k+1)
                    });
    }
    return std::unique_ptr<Bulk>(factory.createGrid().release());
}

//! a helix and straight vessels in [0,1]^3
std::unique_ptr<Network> makeNetwork()
{
    Dune::GridFactory<Network> factory;
    unsigned int numVertices = 0;
    const auto addPolyline = [&](const std::vector<Point>& points)
    {
        for (std::size_t i = 0; i < points.size(); ++i)
        {
            factory.insertVertex(points[i]);
            if (i > 0)
                factory.insertElement(Dune::GeometryTypes::line, {numVertices - 1, numVertices});
            ++numVertices;
        }
    };

    std::vector<Point> helix;
    for (int i = 0; i <= 200; ++i)
    {
        const double t = i/200.0;
        const double phi = 0.37 + 2.0*M_PI*3.3*t;
        helix.push_back({0.5 + 0.35*std::cos(phi), 0.5 + 0.35*std::sin(phi), 0.06 + 0.88*t});
    }
    addPolyline(helix);

    for (int v = 0; v < 6; ++v)
    {
        std::vector<Point> straight;
        const Point a{0.1 + 0.15*v, 0.05, 0.1}, b{0.9 - 0.1*v, 0.95, 0.9};
        for (int i = 0; i <= 60; ++i)
        {
            auto x = a;
            x.axpy(i/60.0, b - a);
            straight.push_back(x);
        }
        addPolyline(straight);
    }
    return factory.createGrid();
}

//! the weight of each element of this process, for ALUGrid's partitioner
struct ElementWeights
{
    const Bulk& grid;
    const std::vector<int>& weights;
    int operator()(const Bulk::Codim<0>::Entity& element) const
    { return weights[grid.levelIndexSet(0).index(element)]; }
};

//! the rank each element of this process goes to, for ALUGrid's repartition
struct ElementDestinations
{
    const Bulk& grid;
    const std::vector<int>& destinations;
    int operator()(const Bulk::Codim<0>::Entity& element) const
    { return destinations[grid.levelIndexSet(0).index(element)]; }
    bool importRanks(std::set<int>&) const
    { return false; }
};

//! the face-neighbor graph of the elements of the coarsest level
std::vector<std::vector<int>> elementGraph(const Bulk& grid)
{
    const auto gridView = grid.levelGridView(0);
    std::vector<std::vector<int>> graph(gridView.size(0));
    for (const auto& element : elements(gridView))
        for (const auto& intersection : intersections(gridView, element))
            if (intersection.neighbor())
                graph[gridView.indexSet().index(element)].push_back(gridView.indexSet().index(intersection.outside()));
    return graph;
}

} // end namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    Parameters::init(argc, argv);

    const auto method = getParam<std::string>("LoadBalance.Method");
    const auto pieceWeight = getParam<double>("LoadBalance.PieceWeight");
    const auto maxImbalance = getParam<double>("LoadBalance.MaxImbalance");

    auto bulk = makeBulk(getParam<int>("Grid.Cells"), mpiHelper.rank());
    bulk->globalRefine(getParam<int>("Grid.Refinement", 0));
    const auto network = makeNetwork();
    const auto& comm = bulk->comm();

    // the weights of the undistributed grid
    auto weights = bulkLoadWeights(bulk->leafGridView(), network->leafGridView(), pieceWeight);
    long totalWeight = 0;
    for (const auto w : weights)
        totalWeight += w;
    totalWeight = comm.sum(totalWeight);

    if (method == "ALUGrid")
    {
        ElementWeights elementWeights{*bulk, weights};
        bulk->loadBalance(elementWeights);
    }
    else if (method == "Scotch")
    {
        std::vector<int> destinations;
        if (mpiHelper.rank() == 0)
            destinations = ScotchPartitioner<int>::partition(elementGraph(*bulk), weights, comm.size());
        ElementDestinations elementDestinations{*bulk, destinations};
        bulk->repartition(elementDestinations);
    }
    else
        bulk->loadBalance();

    // the weights of the distributed grid
    weights = bulkLoadWeights(bulk->leafGridView(), network->leafGridView(), pieceWeight);
    long localWeight = 0;
    for (const auto w : weights)
        localWeight += w;
    const auto distributedWeight = comm.sum(localWeight);
    const auto maxWeight = comm.max(localWeight);
    const double imbalance = maxWeight/(double(totalWeight)/comm.size());

    if (mpiHelper.rank() == 0)
        std::cout << method << " on " << comm.size() << " processes: total weight " << totalWeight
                  << ", after the distribution " << distributedWeight << ", largest on a process " << maxWeight
                  << " (" << imbalance << " times the mean)" << std::endl;

    int failures = 0;
    if (distributedWeight != totalWeight)
    {
        if (mpiHelper.rank() == 0)
            std::cerr << "The weights of the distributed grid do not add up to the weights of the whole grid" << std::endl;
        ++failures;
    }
    if (imbalance > maxImbalance)
    {
        if (mpiHelper.rank() == 0)
            std::cerr << "The largest weight on a process exceeds the mean by more than a factor " << maxImbalance << std::endl;
        ++failures;
    }

    return failures == 0 ? 0 : 1;
}
