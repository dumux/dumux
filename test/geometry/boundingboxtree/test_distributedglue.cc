//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Test the distributed (MPI-parallel) intersection of two bounding box trees.
//
// We intersect two independently partitioned structured grids of the same unit
// cube. The total volume of all pairwise cell intersections must equal the
// volume of the unit cube (1.0), independent of the number of processes. This is
// a strong, partition-independent correctness check for the distributed glue.
//
#include <config.h>

#include <iostream>
#include <vector>
#include <array>
#include <bitset>
#include <cmath>
#include <memory>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/distributedboundingboxtree.hh>
#include <dumux/geometry/distributedintersectingentities.hh>
#include <dumux/geometry/distributedintersectionentityset.hh>
#include <dumux/discretization/basicgridgeometry.hh>
#include <dumux/multidomain/distributedglue.hh>

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
    using EntitySet = GridViewGeometricEntitySet<GridView, 0>;
    using DistTree = DistributedBoundingBoxTree<EntitySet>;

    // two independently partitioned grids of the same unit cube
    Dune::FieldVector<double, dimWorld> upperRight(1.0);
    std::array<int, dimWorld> cellsA; cellsA.fill(8);
    std::array<int, dimWorld> cellsB; cellsB.fill(4);
    Grid gridA(upperRight, cellsA, std::bitset<dimWorld>(0ULL), /*overlap=*/1);
    Grid gridB(upperRight, cellsB, std::bitset<dimWorld>(0ULL), /*overlap=*/1);

    DistTree treeA(std::make_shared<EntitySet>(gridA.leafGridView()));
    DistTree treeB(std::make_shared<EntitySet>(gridB.leafGridView()));

    // compute the distributed intersections (owned by this process)
    const auto rawIntersections = intersectingEntities(treeA, treeB);

    int failed = 0;

    // every returned intersection must be owned by this process as the domain entity
    for (const auto& is : rawIntersections)
        if (is.domainRank != rank)
        {
            std::cerr << "[rank " << rank << "] FAILED: intersection domain not owned by this rank\n";
            ++failed;
        }

    // sum the volume of all intersection simplices over all processes
    double localVolume = 0.0;
    for (const auto& is : rawIntersections)
    {
        using Geometry = Dune::MultiLinearGeometry<double, dimWorld, dimWorld>;
        Geometry geo(Dune::GeometryTypes::simplex(dimWorld), is.corners);
        localVolume += geo.volume();
    }
    const double globalVolume = comm.sum(localVolume);
    const std::size_t globalCount = comm.sum(rawIntersections.size());

    const double expectedVolume = 1.0; // volume of the unit cube
    if (std::abs(globalVolume - expectedVolume) > 1e-8)
    {
        if (rank == 0)
            std::cerr << "FAILED: total intersection volume " << globalVolume
                      << " != expected " << expectedVolume << "\n";
        ++failed;
    }

    // exercise the DistributedIntersectionEntitySet wrapper and check consistency
    {
        using EntitySetPtr = std::shared_ptr<EntitySet>;
        EntitySetPtr setA = std::make_shared<EntitySet>(gridA.leafGridView());
        EntitySetPtr setB = std::make_shared<EntitySet>(gridB.leafGridView());
        DistributedIntersectionEntitySet<EntitySet, EntitySet> glue;
        glue.build(setA, setB);

        if (glue.size() != rawIntersections.size())
        {
            std::cerr << "[rank " << rank << "] FAILED: wrapper size " << glue.size()
                      << " != free-function size " << rawIntersections.size() << "\n";
            ++failed;
        }

        double wrapperVolume = 0.0;
        for (const auto& is : intersections(glue))
        {
            wrapperVolume += is.geometry().volume();
            if (is.domainRank() != rank)
            {
                std::cerr << "[rank " << rank << "] FAILED: wrapper intersection domain not owned\n";
                ++failed;
            }
            // the domain entity must be accessible and interior on this rank
            if (is.domainEntity().partitionType() != Dune::InteriorEntity)
            {
                std::cerr << "[rank " << rank << "] FAILED: wrapper domain entity not interior\n";
                ++failed;
            }
        }
        const double globalWrapperVolume = comm.sum(wrapperVolume);
        if (std::abs(globalWrapperVolume - expectedVolume) > 1e-8)
        {
            if (rank == 0)
                std::cerr << "FAILED: wrapper total intersection volume " << globalWrapperVolume
                          << " != expected " << expectedVolume << "\n";
            ++failed;
        }
    }

    // exercise the makeDistributedGlue convenience using grid geometries
    {
        using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
        using GridGeometry = BasicGridGeometry<GridView, Mapper, Mapper>;
        GridGeometry ggA(gridA.leafGridView());
        GridGeometry ggB(gridB.leafGridView());

        const auto glue = makeDistributedGlue(ggA, ggB);
        double glueVolume = 0.0;
        for (const auto& is : intersections(glue))
            glueVolume += is.geometry().volume();
        const double globalGlueVolume = comm.sum(glueVolume);
        if (std::abs(globalGlueVolume - expectedVolume) > 1e-8)
        {
            if (rank == 0)
                std::cerr << "FAILED: makeDistributedGlue total volume " << globalGlueVolume
                          << " != expected " << expectedVolume << "\n";
            ++failed;
        }
    }

    // exercise the domain interior+overlap partition policy: on top of the interior
    // intersections this additionally reports, on each process, the intersections of the
    // domain entities it overlaps (the intended redundancy for an overlapping assembly).
    // Two checks: (a) only interior/overlap domain entities may be reported, and (b)
    // filtering the result back to interior domain entities must recover exactly the
    // default tessellation, i.e. its total volume is again 1.0 (the overlap entries are
    // pure redundancy and must not be counted).
    {
        const auto overlapInts = intersectingEntities(treeA, treeB,
            Dune::Partitions::interior + Dune::Partitions::overlap);

        double interiorVolume = 0.0;
        for (const auto& is : overlapInts)
        {
            const auto pt = treeA.entitySet().entity(is.domainIndex).partitionType();
            if (pt != Dune::InteriorEntity && pt != Dune::OverlapEntity)
            {
                std::cerr << "[rank " << rank << "] FAILED: overlap-policy domain entity "
                             "neither interior nor overlap\n";
                ++failed;
            }
            if (pt == Dune::InteriorEntity)
            {
                using Geometry = Dune::MultiLinearGeometry<double, dimWorld, dimWorld>;
                Geometry geo(Dune::GeometryTypes::simplex(dimWorld), is.corners);
                interiorVolume += geo.volume();
            }
        }

        const double globalInteriorVolume = comm.sum(interiorVolume);
        if (std::abs(globalInteriorVolume - expectedVolume) > 1e-8)
        {
            if (rank == 0)
                std::cerr << "FAILED: interior+overlap policy interior-domain volume "
                          << globalInteriorVolume << " != expected " << expectedVolume << "\n";
            ++failed;
        }

        // the interior+overlap result is a superset of the interior-only result
        const std::size_t globalOverlapCount = comm.sum(overlapInts.size());
        if (globalOverlapCount < globalCount)
        {
            if (rank == 0)
                std::cerr << "FAILED: interior+overlap count " << globalOverlapCount
                          << " < interior-only count " << globalCount << "\n";
            ++failed;
        }
        if (rank == 0)
            std::cout << "Distributed glue interior+overlap policy (dim=" << dimWorld
                      << ", ranks=" << numProc << "): " << globalOverlapCount
                      << " intersections (interior-only: " << globalCount << ")\n";
    }

    const int globalFailed = comm.sum(failed);
    if (rank == 0)
    {
        std::cout << "Distributed glue (dim=" << dimWorld << ", ranks=" << numProc << "): "
                  << globalCount << " intersections, total volume = " << globalVolume << "\n";
        if (globalFailed == 0)
            std::cout << "All distributed glue tests passed." << std::endl;
        else
            std::cerr << globalFailed << " distributed glue test(s) failed." << std::endl;
    }

    return globalFailed != 0 ? 1 : 0;
}
