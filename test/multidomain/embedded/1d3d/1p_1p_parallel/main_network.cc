// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief MPI-parallel 1d-3d embedded coupling on a realistic DGF network (lupine root system).
 *
 * Demonstrates the DGF-network capability of the overlapping (Strategy A) decomposition: the full
 * network is read replicated, its per-element radius is computed from the DGF surface parameter,
 * the network is distributed by the bulk YaspGrid partition, and the radius is migrated onto the
 * local elements by global id and injected into the spatial parameters. The distributed coupled
 * solution is compared, owned dof by owned dof, against a replicated sequential reference; the
 * migrated radii are also checked against the reference by global id.
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
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/yaspgrid/partitioning.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/dgffoam.hh>
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

#include "properties_network.hh"

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
    using LowDimSpatialParams = GetPropType<LowDimTypeTag, Properties::SpatialParams>;
    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    using LATraits = LinearAlgebraTraitsFromAssembler<Assembler>;
    using DirectSolver = DirectSolverMumps<SeqLinearSolverTraits, LATraits>;

    using GlobalPosition = Dune::FieldVector<double, 3>;
    const auto bulkLower = getParam<GlobalPosition>("Tissue.Grid.LowerLeft");
    const auto bulkUpper = getParam<GlobalPosition>("Tissue.Grid.UpperRight");
    const auto bulkCells = getParam<std::array<int, 3>>("Tissue.Grid.Cells");
    const int overlap = getParam<int>("Tissue.Grid.Overlap", 1);
    const auto networkFile = getParam<std::string>("Vessel.Grid.File");
    const int surfaceIdx = getParam<int>("Vessel.Grid.SurfaceParamIndex", 2);

    std::array<double, 3> h{};
    for (int d = 0; d < 3; ++d)
        h[d] = (bulkUpper[d] - bulkLower[d]) / bulkCells[d];

    const auto bulkCanonicalIndex = [&](const GlobalPosition& center) -> std::size_t
    {
        std::array<int, 3> ijk{};
        for (int d = 0; d < 3; ++d)
            ijk[d] = std::clamp(int((center[d] - bulkLower[d]) / h[d]), 0, bulkCells[d] - 1);
        return (std::size_t(ijk[2]) * bulkCells[1] + ijk[1]) * bulkCells[0] + ijk[0];
    };

    // per-element radius from the DGF surface parameter: r = surface / length / (2 pi)
    const auto radiiFromGridView = [&](const auto& gridView, const auto& parameters) -> std::vector<double>
    {
        std::vector<double> radii(gridView.size(0));
        for (const auto& e : elements(gridView))
        {
            const auto eIdx = gridView.indexSet().index(e);
            const double surface = parameters(e)[surfaceIdx];
            const double length = e.geometry().volume();
            radii[eIdx] = surface / length / (2.0 * M_PI);
        }
        return radii;
    };

    // the coupled stationary 1p solve (linear -> single assemble + solve)
    const auto solveCoupled = [&](std::shared_ptr<BulkGridGeometry> bulkGG,
                                  std::shared_ptr<LowDimGridGeometry> lowDimGG,
                                  std::shared_ptr<LowDimSpatialParams> lowDimSpatialParams) -> SolutionVector
    {
        auto couplingManager = std::make_shared<CouplingManager>(bulkGG, lowDimGG);
        auto bulkProblem = std::make_shared<BulkProblem>(bulkGG, couplingManager);
        auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimGG, lowDimSpatialParams, couplingManager);

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
    // 1) reference: full replicated network + bulk on a self communicator
    // ------------------------------------------------------------------
    const auto selfComm = typename BulkGrid::Communication(mpiHelper.getLocalCommunicator());
    BulkGrid bulkRefGrid(bulkLower, bulkUpper, bulkCells, std::bitset<3>(0), overlap, selfComm);

    // read the full network replicated (self communicator) and its DGF parameters
    Dune::GridPtr<LowDimGrid> networkPtr(networkFile, mpiHelper.getLocalCommunicator());
    auto& fullNetwork = *networkPtr;
    const auto dgfParams = [&](const auto& e) -> const std::vector<double>& { return networkPtr.parameters(e); };
    const auto fullRadii = radiiFromGridView(fullNetwork.leafGridView(), dgfParams);

    auto bulkRefGG = std::make_shared<BulkGridGeometry>(bulkRefGrid.leafGridView());
    auto lowRefGG = std::make_shared<LowDimGridGeometry>(fullNetwork.leafGridView());
    auto lowRefSpatialParams = std::make_shared<LowDimSpatialParams>(lowRefGG, fullRadii);
    const auto solRef = solveCoupled(bulkRefGG, lowRefGG, lowRefSpatialParams);

    std::vector<double> bulkRefByCanon(std::size_t(bulkCells[0]) * bulkCells[1] * bulkCells[2], 0.0);
    {
        const auto& gv = bulkRefGrid.leafGridView();
        const auto& is = gv.indexSet();
        for (const auto& e : elements(gv))
            bulkRefByCanon[bulkCanonicalIndex(e.geometry().center())] = solRef[bulkIdx][is.index(e)][0];
    }
    const auto& lowRef = solRef[lowDimIdx]; // indexed by full-network element index == global id

    // ------------------------------------------------------------------
    // 2) parallel: bulk split along z; network distributed by (nested in) the bulk partition
    // ------------------------------------------------------------------
    const auto worldComm = typename BulkGrid::Communication(mpiHelper.getCommunicator());
    const auto partitioning = getParam<std::array<int, 3>>(
        "Tissue.Grid.Partitioning", std::array<int, 3>{1, 1, worldComm.size()});
    if (partitioning[0] * partitioning[1] * partitioning[2] != worldComm.size())
        DUNE_THROW(Dune::InvalidStateException, "Tissue.Grid.Partitioning product must equal the number of ranks");
    Dune::Yasp::FixedSizePartitioning<3> partitioner(partitioning);
    BulkGrid bulkParGrid(bulkLower, bulkUpper, bulkCells, std::bitset<3>(0), overlap, worldComm, &partitioner);

    Dune::Communication<Dune::MPIHelper::MPICommunicator> comm(mpiHelper.getCommunicator());
    auto localNetwork = Dune::FoamGridParallel::distributeFromSpatialPartition(
        fullNetwork, bulkParGrid.leafGridView(), comm);
    const auto par = localNetwork->parallelData();

    // migrate the per-element radius onto the local network elements by global id
    std::vector<double> localRadii(localNetwork->leafGridView().size(0));
    {
        const auto& gv = localNetwork->leafGridView();
        const auto& is = gv.indexSet();
        for (const auto& e : elements(gv))
            localRadii[is.index(e)] = fullRadii[par->elementGlobalId(is.index(e))];
    }

    auto bulkParGG = std::make_shared<BulkGridGeometry>(bulkParGrid.leafGridView());
    auto lowParGG = std::make_shared<LowDimGridGeometry>(localNetwork->leafGridView());
    auto lowParSpatialParams = std::make_shared<LowDimSpatialParams>(lowParGG, localRadii);
    const auto solPar = solveCoupled(bulkParGG, lowParGG, lowParSpatialParams);

    // ------------------------------------------------------------------
    // 3) compare owned dofs (and migrated radii) against the reference by global id
    // ------------------------------------------------------------------
    double maxDiff = 0.0, maxRef = 0.0, maxRadiusDiff = 0.0;
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
        for (const auto& e : elements(gv, Dune::Partitions::interior))
        {
            const auto gid = par->elementGlobalId(is.index(e));
            maxDiff = std::max(maxDiff, std::abs(solPar[lowDimIdx][is.index(e)][0] - lowRef[gid][0]));
            maxRef = std::max(maxRef, std::abs(lowRef[gid][0]));
            maxRadiusDiff = std::max(maxRadiusDiff, std::abs(localRadii[is.index(e)] - fullRadii[gid]));
            ++numChecked;
        }
    }

    maxDiff = comm.max(maxDiff);
    maxRef = comm.max(maxRef);
    maxRadiusDiff = comm.max(maxRadiusDiff);
    const std::size_t totalChecked = comm.sum(numChecked);
    const double relErr = maxDiff / (maxRef > 0.0 ? maxRef : 1.0);
    const bool ok = (relErr < 1e-7) && (maxRadiusDiff == 0.0);

    if (comm.rank() == 0)
    {
        std::cout << "1d3d DGF network (cc): compared " << totalChecked << " interior dofs over "
                  << comm.size() << " rank(s); max |parallel - sequential| = " << maxDiff
                  << " (relative " << relErr << "), max radius mismatch = " << maxRadiusDiff << "\n";
        std::cout << (ok ? "PASS: parallel DGF-network solution and migrated radii match sequential\n"
                         : "FAIL: parallel DGF-network solution/radii differ from sequential\n");
    }

    return ok ? 0 : 1;
}
