// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Stationary 1p flow on a distributed network: the parallel solution equals the sequential one
 *
 * Every rank solves on the whole network as the reference, then on a copy distributed either by a
 * YaspGrid around the network (as a bulk grid would distribute it) or by coordinate bisection, and
 * compares all local dofs with the reference through the global ids.
 */
#include <config.h>

#include <array>
#include <bitset>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
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
#include "bulkpartition.hh"

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

    // refine elements longer than Grid.MaxElementLength until there are none (collective)
    const auto maxElementLength = getParam<double>("Grid.MaxElementLength", 0.0);
    auto refine = [&](Grid& grid)
    {
        if (maxElementLength <= 0.0)
            return;
        while (true)
        {
            int marked = 0;
            for (const auto& element : elements(grid.leafGridView()))
                if (element.geometry().volume() > maxElementLength)
                {
                    grid.mark(1, element);
                    marked = 1;
                }
            if (!grid.comm().max(marked))
                return;
            grid.preAdapt();
            grid.adapt();
            grid.postAdapt();
        }
    };

    // position of a degree of freedom, rounded so that it identifies the dof in both grids
    auto dofKey = [](const auto& entity)
    {
        const auto x = entity.geometry().center();
        return std::array<long long, 3>{std::llround(x[0]*1e9), std::llround(x[1]*1e9), std::llround(x[2]*1e9)};
    };

    // --- 1) reference: full grid, sequential (self communicator), on every rank ---
    Dune::GridPtr<Grid> fullPtr("network.dgf", mpiHelper.getLocalCommunicator());
    Grid& fullGrid = *fullPtr; // a plain FoamGrid carries a self communicator -> sequential solve
    refine(fullGrid);
    const auto xFull = solve(fullGrid);
    std::map<std::array<long long, 3>, double> reference;
    for (const auto& entity : entities(fullGrid.leafGridView(), Dune::Codim<dofCodim>{}))
        reference[dofKey(entity)] = xFull[fullGrid.leafGridView().indexSet().index(entity)][0];

    // --- 2) parallel: distribute a copy of the network over the ranks and solve in parallel ---
    Dune::Communication<Dune::MPIHelper::MPICommunicator> comm(mpiHelper.getCommunicator());
    Dune::GridPtr<Grid> localPtr("network.dgf", mpiHelper.getLocalCommunicator());
    Grid& localGrid = *localPtr;
    localGrid.setCommunicator(mpiHelper.getCommunicator());
    const auto overlap = getParam<int>("Grid.Overlap", 1);
    if (getParam<bool>("Grid.PartitionByBulk", true))
        localGrid.loadBalance(partitionByBoxAround(localGrid, getParam<int>("Grid.BulkCells", 8)), overlap);
    else
        localGrid.loadBalance(Dune::FoamGridParallel::recursiveCoordinateBisection(localGrid.leafGridView(), comm.size()), overlap);

    refine(localGrid);
    const auto localGV = localGrid.leafGridView();
    const auto xLocal = solve(localGrid);

    // --- 3) compare every local dof (interior/border/overlap copies) with the reference at the
    //         same position. After a correct parallel solve all copies hold their owner's value. ---
    double maxDiff = 0.0, maxRef = 0.0;
    int numChecked = 0, numMissing = 0;
    for (const auto& entity : entities(localGV, Dune::Codim<dofCodim>{}))
    {
        const auto it = reference.find(dofKey(entity));
        if (it == reference.end())
        {
            ++numMissing;
            continue;
        }
        const auto localIndex = localGV.indexSet().index(entity);
        maxDiff = std::max(maxDiff, std::abs(xLocal[localIndex][0] - it->second));
        maxRef = std::max(maxRef, std::abs(it->second));
        ++numChecked;
    }
    numMissing = comm.sum(numMissing);
    maxDiff = comm.max(maxDiff);
    maxRef = comm.max(maxRef);
    const int totalChecked = comm.sum(numChecked);
    const double relErr = maxDiff / (maxRef > 0.0 ? maxRef : 1.0);

    if (comm.rank() == 0)
    {
        std::cout << "1p network (" << (dofCodim == 0 ? "cc" : "box") << "): compared "
                  << totalChecked << " local dofs over " << comm.size()
                  << " rank(s), " << fullGrid.leafGridView().size(dofCodim) << " dofs in the whole grid"
                  << "; max |parallel - sequential| = " << maxDiff
                  << " (relative " << relErr << ")\n";
    }

    const bool ok = (relErr < 1e-7) && numMissing == 0 && totalChecked > 0;
    if (comm.rank() == 0 && numMissing > 0)
        std::cout << numMissing << " local dofs have no counterpart in the sequential grid\n";
    if (comm.rank() == 0)
        std::cout << (ok ? "PASS: parallel 1p solution matches sequential\n"
                         : "FAIL: parallel 1p solution differs from sequential\n");
    return ok ? 0 : 1;
}
