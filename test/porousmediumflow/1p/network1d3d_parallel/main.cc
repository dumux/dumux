// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Validate the distributed (overlapping) parallel FoamGrid: solve a static 1p cc-tpfa
 *        problem on a 1d network and check that the parallel solution matches the sequential one.
 *
 * Every rank loads the full network (replicated, self-communicator) and solves it sequentially as
 * the reference. The network is then distributed by the spatial decomposition of a coarse YaspGrid
 * and solved in parallel. Each rank compares its owned (interior) dofs against the reference using
 * the grid's global ids.
 */
#include <config.h>

#include <array>
#include <bitset>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

#include <set>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include "properties.hh"

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>

#ifndef TYPETAG
#define TYPETAG OnePNetworkCCTpfa
#endif

#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/dgffoam.hh>
#include <dune/foamgrid/parallel/distribute.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TYPETAG;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    // dof codimension: 0 for cell-centered, dim (vertices) for box
    static constexpr int dofCodim =
        (GridGeometry::discMethod == DiscretizationMethods::box) ? int(Grid::dimension) : 0;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;

    // stationary 1p solve on a given grid; returns the solution indexed by leaf element index
    auto solve = [](const Grid& grid) -> SolutionVector
    {
        const auto gridView = grid.leafGridView();
        auto gridGeometry = std::make_shared<GridGeometry>(gridView);
        auto problem = std::make_shared<Problem>(gridGeometry);
        SolutionVector x(gridGeometry->numDofs());
        problem->applyInitialSolution(x);
        auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
        gridVariables->init(x);
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
        auto linearSolver = std::make_shared<LinearSolver>(gridView, gridGeometry->dofMapper());
        NewtonSolver<Assembler, LinearSolver> newton(assembler, linearSolver);
        newton.solve(x);
        return x;
    };

    // --- 1) reference: full grid, sequential (self communicator), on every rank ---
    Dune::GridPtr<Grid> fullPtr("network.dgf", mpiHelper.getLocalCommunicator());
    Grid& fullGrid = *fullPtr; // a plain FoamGrid carries a self communicator -> sequential solve
    const auto xFull = solve(fullGrid);

    // --- 2) parallel: distribute by a coarse YaspGrid, solve in parallel ---
    Dune::Communication<Dune::MPIHelper::MPICommunicator> comm(mpiHelper.getCommunicator());

    Dune::FieldVector<double, 3> lo(std::numeric_limits<double>::max());
    Dune::FieldVector<double, 3> hi(std::numeric_limits<double>::lowest());
    for (const auto& v : vertices(fullGrid.leafGridView()))
    {
        const auto x = v.geometry().center();
        for (int d = 0; d < 3; ++d) { lo[d] = std::min(lo[d], x[d]); hi[d] = std::max(hi[d], x[d]); }
    }
    for (int d = 0; d < 3; ++d) { const double m = 0.05 * (hi[d] - lo[d]) + 1e-12; lo[d] -= m; hi[d] += m; }

    using Coarse = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>;
    const int coarseCells = getParam<int>("Grid.CoarseCells", 8);
    const int coarseOverlap = getParam<int>("Grid.CoarseOverlap", 1);
    Coarse coarse(lo, hi, std::array<int, 3>{coarseCells, coarseCells, coarseCells},
                  std::bitset<3>(0), coarseOverlap, comm);

    auto local = Dune::FoamGridParallel::distributeFromSpatialPartition(fullGrid, coarse.leafGridView(), comm);
    const auto localGV = local->leafGridView();
    const auto par = local->parallelData();

    // --- diagnostic: are all interior elements' face-neighbours present locally? ---
    {
        const auto& fGV = fullGrid.leafGridView();
        const auto& fIS = fGV.indexSet();
        std::vector<std::vector<std::size_t>> facetElems(fIS.size(1));
        for (const auto& e : elements(fGV))
        {
            const auto refE = referenceElement(e.geometry());
            for (int i = 0; i < refE.size(1); ++i)
                facetElems[fIS.subIndex(e, i, 1)].push_back(fIS.index(e));
        }
        std::vector<std::vector<std::size_t>> adj(fIS.size(0));
        for (const auto& sh : facetElems)
            for (std::size_t a = 0; a < sh.size(); ++a)
                for (std::size_t b = a + 1; b < sh.size(); ++b)
                { adj[sh[a]].push_back(sh[b]); adj[sh[b]].push_back(sh[a]); }

        std::set<std::size_t> localGids;
        for (const auto& e : elements(localGV))
            localGids.insert(par->elementGlobalId(localGV.indexSet().index(e)));

        int missing = 0;
        for (const auto& e : elements(localGV))
        {
            if (e.partitionType() != Dune::InteriorEntity) continue;
            const auto gid = par->elementGlobalId(localGV.indexSet().index(e));
            for (const auto n : adj[gid]) if (!localGids.count(n)) ++missing;
        }
        const int totalMissing = comm.sum(missing);
        if (comm.rank() == 0)
            std::cout << "stencil check (overlap=" << coarseOverlap << "): " << totalMissing
                      << " missing interior face-neighbours\n";
    }

    const auto xLocal = solve(*local);

    // --- 3) compare every local dof (interior/border/overlap copies) against the reference,
    //         by global id. After a correct parallel solve all copies hold their owner's value. ---
    double maxDiff = 0.0, maxRef = 0.0;
    int numChecked = 0;
    for (const auto& entity : entities(localGV, Dune::Codim<dofCodim>{}))
    {
        const auto localIndex = localGV.indexSet().index(entity);
        std::size_t globalId;
        if constexpr (dofCodim == 0)
            globalId = par->elementGlobalId(localIndex);
        else
            globalId = par->vertexGlobalId(localIndex);
        maxDiff = std::max(maxDiff, std::abs(xLocal[localIndex][0] - xFull[globalId][0]));
        maxRef = std::max(maxRef, std::abs(xFull[globalId][0]));
        ++numChecked;
    }
    maxDiff = comm.max(maxDiff);
    maxRef = comm.max(maxRef);
    const int totalChecked = comm.sum(numChecked);
    const double relErr = maxDiff / (maxRef > 0.0 ? maxRef : 1.0);

    if (comm.rank() == 0)
    {
        std::cout << "1p network (" << (dofCodim == 0 ? "cc" : "box") << "): compared "
                  << totalChecked << " local dofs over " << comm.size()
                  << " rank(s); max |parallel - sequential| = " << maxDiff
                  << " (relative " << relErr << ")\n";
    }

    const bool ok = (relErr < 1e-7);
    if (comm.rank() == 0)
        std::cout << (ok ? "PASS: parallel 1p solution matches sequential\n"
                         : "FAIL: parallel 1p solution differs from sequential\n");
    return ok ? 0 : 1;
}
