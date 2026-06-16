//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>
#include <vector>
#include <array>
#include <bitset>
#include <string>
#include <memory>
#include <algorithm>
#include <ranges>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/distributedboundingboxtree.hh>
#include <dumux/geometry/distributedintersectingentities.hh>

#ifndef WORLD_DIMENSION
#define WORLD_DIMENSION 3
#endif

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dune::MPIHelper::instance(argc, argv);
    const auto& comm = Dune::MPIHelper::getCommunication();
    const int rank = comm.rank();
    const int numProc = comm.size();

    static constexpr int dimWorld = WORLD_DIMENSION;
    using Grid = Dune::YaspGrid<dimWorld>;
    using GridView = typename Grid::LeafGridView;
    using GlobalPosition = Dune::FieldVector<double, dimWorld>;
    using EntitySet = GridViewGeometricEntitySet<GridView, 0>;
    using DistTree = DistributedBoundingBoxTree<EntitySet>;

    // a structured grid over the unit cube partitioned across all processes
    const int cellsPerDim = 8;
    Dune::FieldVector<double, dimWorld> upperRight(1.0);
    std::array<int, dimWorld> cells; cells.fill(cellsPerDim);
    Grid grid(upperRight, cells, std::bitset<dimWorld>(0ULL), /*overlap=*/1);
    const auto gridView = grid.leafGridView();

    // build the distributed bounding box tree
    DistTree tree(std::make_shared<EntitySet>(gridView));

    int failed = 0;
    const auto check = [&](bool ok, const std::string& msg)
    {
        if (!ok) { std::cerr << "[rank " << rank << "] FAILED: " << msg << std::endl; ++failed; }
    };

    check(tree.comm().size() == numProc, "comm size mismatch");
    check(tree.rank() == rank, "rank mismatch");

    // sample points at cell centers: each is contained in exactly one cell globally
    const double h = 1.0/cellsPerDim;
    std::vector<GlobalPosition> samplePoints;
    {
        // a handful of interior cell centers along the diagonal and axes
        for (int i = 0; i < cellsPerDim; ++i)
        {
            GlobalPosition p(0.0);
            for (int d = 0; d < dimWorld; ++d)
                p[d] = (i + 0.5)*h;
            samplePoints.push_back(p);
        }
        // mixed-coordinate centers
        GlobalPosition p(0.0);
        for (int d = 0; d < dimWorld; ++d)
            p[d] = ((d % cellsPerDim) + 0.5)*h;
        samplePoints.push_back(p);
    }

    for (const auto& p : samplePoints)
    {
        // process collisions are replicated and must be identical on all ranks:
        // query the process tree and map the resulting leaves to their owning ranks
        // (a lazy range over the leaf indices is enough; we only search it below)
        const auto leaves = intersectingEntities(p, tree.processTree());
        const auto candidates = std::views::transform(leaves,
            [&](std::size_t leafIdx){ return tree.processForLeaf(leafIdx); });

        // local owned intersections (rank, localIdx)
        const auto local = intersectingEntities(p, tree);

        // every local hit's owning rank must appear among the candidate processes
        for (const auto& [r, idx] : local)
            check(std::ranges::find(candidates, r) != candidates.end(),
                  "owned hit on a process not reported by the process tree");

        // exactly one owned cell globally contains a cell-center point
        const std::size_t localCount = local.size();
        const auto globalCount = comm.sum(localCount);
        check(globalCount == 1, "cell-center point not contained in exactly one owned cell globally");
    }

    // a point well outside the domain should not be contained anywhere
    {
        GlobalPosition outside(10.0);
        const auto local = intersectingEntities(outside, tree);
        const auto globalCount = comm.sum(local.size());
        check(globalCount == 0, "point outside domain reported as contained");
    }

    const int globalFailed = comm.sum(failed);
    if (rank == 0)
    {
        if (globalFailed == 0)
            std::cout << "All distributed bounding box tree tests passed (dim=" << dimWorld
                      << ", ranks=" << numProc << ")." << std::endl;
        else
            std::cerr << globalFailed << " distributed bounding box tree test(s) failed." << std::endl;
    }

    return globalFailed != 0 ? 1 : 0;
}
