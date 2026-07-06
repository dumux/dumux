// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief MPI-parallel test of the 1d-3d embedded coupling (overlapping decomposition).
 *
 * Strategy A: distribute both grids from one spatial partition (the bulk YaspGrid drives the
 * network FoamGrid partition), then run the existing serial coupling manager per-rank on the
 * local grid views and solve the monolithic system with the parallel multidomain MUMPS backend.
 *
 * Every rank also builds the full, replicated problem (bulk YaspGrid + full FoamGrid on a self
 * communicator) and solves it sequentially as the reference. The distributed parallel solution
 * is then compared, owned (interior) dof by owned dof, against the reference: bulk cells matched
 * by their structured (canonical) index, network elements by the FoamGrid global id. A correct
 * parallel solve reproduces the sequential one to machine precision.
 */
#include <config.h>

#include <array>
#include <bitset>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/yaspgrid/partitioning.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/parallel/distribute.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/mumpssolver.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    Parameters::init(argc, argv);

    using BulkTypeTag = Properties::TTag::BULKTYPETAG;
    using LowDimTypeTag = Properties::TTag::LOWDIMTYPETAG;

    using BulkGrid = GetPropType<BulkTypeTag, Properties::Grid>;     // YaspGrid<3>
    using LowDimGrid = GetPropType<LowDimTypeTag, Properties::Grid>; // FoamGrid<1,3>
    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using LowDimGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;

    using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    constexpr auto bulkIdx = Traits::template SubDomain<0>::Index();
    constexpr auto lowDimIdx = Traits::template SubDomain<1>::Index();
    using SolutionVector = typename Traits::SolutionVector;

    using CouplingManager = GetPropType<BulkTypeTag, Properties::CouplingManager>;
    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    using LATraits = LinearAlgebraTraitsFromAssembler<Assembler>;
    using DirectSolver = DirectSolverMumps<SeqLinearSolverTraits, LATraits>;

    static_assert(BulkGridGeometry::discMethod == DiscretizationMethods::cctpfa
                  && LowDimGridGeometry::discMethod == DiscretizationMethods::cctpfa,
                  "This parallel test currently assumes cell-centered (codim-0) dofs in both domains");

    // ------------------------------------------------------------------
    // grid parameters (shared by the reference and the parallel grids)
    // ------------------------------------------------------------------
    using GlobalPosition = Dune::FieldVector<double, 3>;
    const auto bulkLower = getParam<GlobalPosition>("Tissue.Grid.LowerLeft");
    const auto bulkUpper = getParam<GlobalPosition>("Tissue.Grid.UpperRight");
    const auto bulkCells = getParam<std::array<int, 3>>("Tissue.Grid.Cells");
    const int overlap = getParam<int>("Tissue.Grid.Overlap", 1);

    const auto netLower = getParam<GlobalPosition>("Vessel.Grid.LowerLeft");
    const auto netUpper = getParam<GlobalPosition>("Vessel.Grid.UpperRight");
    const int netSegments = getParam<int>("Vessel.Grid.Cells");

    std::array<double, 3> h{};
    for (int d = 0; d < 3; ++d)
        h[d] = (bulkUpper[d] - bulkLower[d]) / bulkCells[d];

    // canonical (structured, partition-independent) index of a bulk cell from its center
    const auto bulkCanonicalIndex = [&](const GlobalPosition& center) -> std::size_t
    {
        std::array<int, 3> ijk{};
        for (int d = 0; d < 3; ++d)
        {
            int i = static_cast<int>((center[d] - bulkLower[d]) / h[d]);
            ijk[d] = std::clamp(i, 0, bulkCells[d] - 1);
        }
        return (std::size_t(ijk[2]) * bulkCells[1] + ijk[1]) * bulkCells[0] + ijk[0];
    };

    // build the full (replicated) network FoamGrid on every rank: a structured 1d interval
    const auto makeFullNetwork = []( const GlobalPosition& lo, const GlobalPosition& hi, int n)
    {
        Dune::GridFactory<LowDimGrid> factory;
        GlobalPosition step = hi; step -= lo; step /= n;
        GlobalPosition pos = lo;
        for (int v = 0; v <= n; ++v, pos += step)
            factory.insertVertex(pos);
        for (int e = 0; e < n; ++e)
            factory.insertElement(Dune::GeometryTypes::line, {unsigned(e), unsigned(e+1)});
        return std::shared_ptr<LowDimGrid>(factory.createGrid());
    };

    // ------------------------------------------------------------------
    // the coupled stationary 1p solve (linear -> single assemble + solve)
    // ------------------------------------------------------------------
    const auto solveCoupled = [&](std::shared_ptr<BulkGridGeometry> bulkGG,
                                  std::shared_ptr<LowDimGridGeometry> lowDimGG) -> SolutionVector
    {
        auto couplingManager = std::make_shared<CouplingManager>(bulkGG, lowDimGG);
        auto bulkProblem = std::make_shared<BulkProblem>(bulkGG, couplingManager);
        auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimGG, couplingManager);

        SolutionVector sol;
        sol[bulkIdx].resize(bulkGG->numDofs());
        sol[lowDimIdx].resize(lowDimGG->numDofs());
        bulkProblem->applyInitialSolution(sol[bulkIdx]);
        lowDimProblem->applyInitialSolution(sol[lowDimIdx]);

        couplingManager->init(bulkProblem, lowDimProblem, sol);
        bulkProblem->computePointSourceMap();
        lowDimProblem->computePointSourceMap();

        auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkGG);
        bulkGridVariables->init(sol[bulkIdx]);
        auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimGG);
        lowDimGridVariables->init(sol[lowDimIdx]);

        auto assembler = std::make_shared<Assembler>(
            std::make_tuple(bulkProblem, lowDimProblem),
            std::make_tuple(bulkGG, lowDimGG),
            std::make_tuple(bulkGridVariables, lowDimGridVariables),
            couplingManager);

        const bool parallel = bulkGG->gridView().comm().size() > 1;
        std::shared_ptr<DirectSolver> linearSolver = parallel
            ? std::make_shared<DirectSolver>(std::make_tuple(bulkGG, lowDimGG))
            : std::make_shared<DirectSolver>();

        couplingManager->updateSolution(sol);
        assembler->assembleJacobianAndResidual(sol);

        SolutionVector deltaSol = sol;
        linearSolver->solve(assembler->jacobian(), deltaSol, assembler->residual());
        sol -= deltaSol;

        couplingManager->updateSolution(sol);
        return sol;
    };

    // ------------------------------------------------------------------
    // 1) reference: full replicated problem on a self communicator
    // ------------------------------------------------------------------
    const auto selfComm = typename BulkGrid::Communication(mpiHelper.getLocalCommunicator());
    BulkGrid bulkRefGrid(bulkLower, bulkUpper, bulkCells, std::bitset<3>(0), overlap, selfComm);
    auto fullNetwork = makeFullNetwork(netLower, netUpper, netSegments);

    auto bulkRefGG = std::make_shared<BulkGridGeometry>(bulkRefGrid.leafGridView());
    auto lowRefGG = std::make_shared<LowDimGridGeometry>(fullNetwork->leafGridView());
    const auto solRef = solveCoupled(bulkRefGG, lowRefGG);

    // reference lookup tables (indexed independently of the parallel partition)
    std::vector<double> bulkRefByCanon(std::size_t(bulkCells[0]) * bulkCells[1] * bulkCells[2], 0.0);
    {
        const auto& gv = bulkRefGrid.leafGridView();
        const auto& is = gv.indexSet();
        for (const auto& e : elements(gv))
            bulkRefByCanon[bulkCanonicalIndex(e.geometry().center())] = solRef[bulkIdx][is.index(e)][0];
    }
    const auto& lowRef = solRef[lowDimIdx]; // indexed by full-network element index == global id

    // ------------------------------------------------------------------
    // 2) parallel: bulk YaspGrid split along z so the network is genuinely distributed,
    //    network FoamGrid partitioned by (nested in) the bulk partition
    // ------------------------------------------------------------------
    const auto worldComm = typename BulkGrid::Communication(mpiHelper.getCommunicator());
    // default: split along z so the z-aligned network is distributed; overridable so that the
    // average-mode test can instead split across the circle plane (x/y) to stress the overlap.
    const auto partitioning = getParam<std::array<int, 3>>(
        "Tissue.Grid.Partitioning", std::array<int, 3>{1, 1, worldComm.size()});
    if (partitioning[0] * partitioning[1] * partitioning[2] != worldComm.size())
        DUNE_THROW(Dune::InvalidStateException, "Tissue.Grid.Partitioning product must equal the number of ranks");
    Dune::Yasp::FixedSizePartitioning<3> partitioner(partitioning);
    BulkGrid bulkParGrid(bulkLower, bulkUpper, bulkCells, std::bitset<3>(0), overlap, worldComm, &partitioner);

    Dune::Communication<Dune::MPIHelper::MPICommunicator> comm(mpiHelper.getCommunicator());
    auto localNetwork = Dune::FoamGridParallel::distributeFromSpatialPartition(
        *fullNetwork, bulkParGrid.leafGridView(), comm);

    auto bulkParGG = std::make_shared<BulkGridGeometry>(bulkParGrid.leafGridView());
    auto lowParGG = std::make_shared<LowDimGridGeometry>(localNetwork->leafGridView());
    const auto solPar = solveCoupled(bulkParGG, lowParGG);

    // ------------------------------------------------------------------
    // 3) compare owned (interior) dofs against the reference by global id
    // ------------------------------------------------------------------
    double maxDiff = 0.0, maxRef = 0.0;
    std::size_t numChecked = 0;

    {
        const auto& gv = bulkParGrid.leafGridView();
        const auto& is = gv.indexSet();
        for (const auto& e : elements(gv, Dune::Partitions::interior))
        {
            const auto ref = bulkRefByCanon[bulkCanonicalIndex(e.geometry().center())];
            maxDiff = std::max(maxDiff, std::abs(solPar[bulkIdx][is.index(e)][0] - ref));
            maxRef = std::max(maxRef, std::abs(ref));
            ++numChecked;
        }
    }
    {
        const auto& gv = localNetwork->leafGridView();
        const auto& is = gv.indexSet();
        const auto par = localNetwork->parallelData();
        for (const auto& e : elements(gv, Dune::Partitions::interior))
        {
            const auto gid = par->elementGlobalId(is.index(e));
            const auto ref = lowRef[gid][0];
            maxDiff = std::max(maxDiff, std::abs(solPar[lowDimIdx][is.index(e)][0] - ref));
            maxRef = std::max(maxRef, std::abs(ref));
            ++numChecked;
        }
    }

    maxDiff = comm.max(maxDiff);
    maxRef = comm.max(maxRef);
    const std::size_t totalChecked = comm.sum(numChecked);
    const double relErr = maxDiff / (maxRef > 0.0 ? maxRef : 1.0);
    const bool ok = relErr < 1e-7;

    if (comm.rank() == 0)
    {
        std::cout << "1d3d embedded (cc): compared " << totalChecked << " interior dofs over "
                  << comm.size() << " rank(s); max |parallel - sequential| = " << maxDiff
                  << " (relative " << relErr << ")\n";
        std::cout << (ok ? "PASS: parallel 1d3d solution matches sequential\n"
                         : "FAIL: parallel 1d3d solution differs from sequential\n");
    }

    return ok ? 0 : 1;
}
