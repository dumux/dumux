// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief MPI-parallel transient root-soil benchmark (1p2c Richards / 1pnc root), realistic DGF
 *        network, overlapping (Strategy A) decomposition solved with the parallel multidomain
 *        MUMPS backend.
 *
 * The full root network is read replicated, its per-element radius is computed from the DGF
 * surface parameter, the network is distributed by the soil YaspGrid partition, and the radius is
 * migrated onto the local elements by global id. A fixed-time-step transient nonlinear solve is
 * run both replicated-sequentially (reference) and distributed-in-parallel; the final-time
 * solution (all components) and the migrated radii are compared owned dof by owned dof.
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
#include <dumux/common/timeloop.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/mumpssolver.hh>
#if defined(USE_RAS_SOLVER)
// Block-diagonal (Restricted Additive Schwarz) parallel solver: exact local MUMPS block solve
// per rank inside a parallel BiCGSTAB, instead of one global parallel direct factorization.
#include <dumux/linear/mumpsrassolver.hh>
#endif

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    Parameters::init(argc, argv);

    using BulkTypeTag = Properties::TTag::Soil;
    using LowDimTypeTag = Properties::TTag::Root;

    using BulkGrid = GetPropType<BulkTypeTag, Properties::Grid>;     // YaspGrid<3>
    using LowDimGrid = GetPropType<LowDimTypeTag, Properties::Grid>; // FoamGrid<1,3>
    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using LowDimGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;

    using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    constexpr auto bulkIdx = Traits::template SubDomain<0>::Index();
    constexpr auto lowDimIdx = Traits::template SubDomain<1>::Index();
    using Scalar = typename Traits::Scalar;
    using SolutionVector = typename Traits::SolutionVector;

    using CouplingManager = GetPropType<BulkTypeTag, Properties::CouplingManager>;
    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    using LowDimSpatialParams = GetPropType<LowDimTypeTag, Properties::SpatialParams>;
    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    using LATraits = LinearAlgebraTraitsFromAssembler<Assembler>;
#if defined(USE_RAS_SOLVER)
    using DirectSolver = MumpsRASSolver<SeqLinearSolverTraits, LATraits>;
#else
    using DirectSolver = DirectSolverMumps<SeqLinearSolverTraits, LATraits>;
#endif
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, DirectSolver, CouplingManager>;

    using GlobalPosition = Dune::FieldVector<double, 3>;
    const auto soilLower = getParam<GlobalPosition>("Soil.Grid.LowerLeft");
    const auto soilUpper = getParam<GlobalPosition>("Soil.Grid.UpperRight");
    const auto soilCells = getParam<std::array<int, 3>>("Soil.Grid.Cells");
    const int overlap = getParam<int>("Soil.Grid.Overlap", 1);
    const auto networkFile = getParam<std::string>("Root.Grid.File");
    const int surfaceIdx = getParam<int>("Root.Grid.SurfaceParamIndex", 2);

    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    std::array<double, 3> h{};
    for (int d = 0; d < 3; ++d)
        h[d] = (soilUpper[d] - soilLower[d]) / soilCells[d];

    const auto soilCanonicalIndex = [&](const GlobalPosition& center) -> std::size_t
    {
        std::array<int, 3> ijk{};
        for (int d = 0; d < 3; ++d)
            ijk[d] = std::clamp(int((center[d] - soilLower[d]) / h[d]), 0, soilCells[d] - 1);
        return (std::size_t(ijk[2]) * soilCells[1] + ijk[1]) * soilCells[0] + ijk[0];
    };

    const auto radiiFromGridView = [&](const auto& gridView, const auto& parameters) -> std::vector<double>
    {
        std::vector<double> radii(gridView.size(0));
        for (const auto& e : elements(gridView))
            radii[gridView.indexSet().index(e)] = parameters(e)[surfaceIdx] / e.geometry().volume() / (2.0 * M_PI);
        return radii;
    };

    // fixed-time-step transient coupled solve to tEnd; returns the final solution
    const auto solveTransient = [&](std::shared_ptr<BulkGridGeometry> bulkGG,
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
        auto oldSol = sol;

        couplingManager->init(bulkProblem, lowDimProblem, sol);
        bulkProblem->computePointSourceMap();
        lowDimProblem->computePointSourceMap();

        auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkGG);
        bulkGridVariables->init(sol[bulkIdx]);
        auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimGG);
        lowDimGridVariables->init(sol[lowDimIdx]);

        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);

        auto assembler = std::make_shared<Assembler>(
            std::make_tuple(bulkProblem, lowDimProblem),
            std::make_tuple(bulkGG, lowDimGG),
            std::make_tuple(bulkGridVariables, lowDimGridVariables),
            couplingManager, timeLoop, oldSol);

        const bool parallel = bulkGG->gridView().comm().size() > 1;
        std::shared_ptr<DirectSolver> linearSolver = parallel
            ? std::make_shared<DirectSolver>(std::make_tuple(bulkGG, lowDimGG))
            : std::make_shared<DirectSolver>();

        NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager, bulkGG->gridView().comm());

        timeLoop->start();
        while (!timeLoop->finished())
        {
            nonLinearSolver.solve(sol, *timeLoop);
            oldSol = sol;
            bulkGridVariables->advanceTimeStep();
            lowDimGridVariables->advanceTimeStep();
            timeLoop->advanceTimeStep();
            timeLoop->setTimeStepSize(dt); // keep the step size fixed (identical steps ref vs parallel)
        }
        return sol;
    };

    // ------------------------------------------------------------------
    // 1) reference: full replicated network + soil on a self communicator
    // ------------------------------------------------------------------
    const auto selfComm = typename BulkGrid::Communication(mpiHelper.getLocalCommunicator());
    BulkGrid soilRefGrid(soilLower, soilUpper, soilCells, std::bitset<3>(0), overlap, selfComm);

    Dune::GridPtr<LowDimGrid> networkPtr(networkFile, mpiHelper.getLocalCommunicator());
    auto& fullNetwork = *networkPtr;
    const auto dgfParams = [&](const auto& e) -> const std::vector<double>& { return networkPtr.parameters(e); };
    const auto fullRadii = radiiFromGridView(fullNetwork.leafGridView(), dgfParams);

    auto soilRefGG = std::make_shared<BulkGridGeometry>(soilRefGrid.leafGridView());
    auto rootRefGG = std::make_shared<LowDimGridGeometry>(fullNetwork.leafGridView());
    auto rootRefSpatialParams = std::make_shared<LowDimSpatialParams>(rootRefGG, fullRadii);
    const auto solRef = solveTransient(soilRefGG, rootRefGG, rootRefSpatialParams);

    constexpr std::size_t bulkBS = std::decay_t<decltype(solRef[bulkIdx][0])>::size();
    constexpr std::size_t lowBS = std::decay_t<decltype(solRef[lowDimIdx][0])>::size();

    std::vector<std::array<double, bulkBS>> soilRefByCanon(std::size_t(soilCells[0]) * soilCells[1] * soilCells[2]);
    {
        const auto& gv = soilRefGrid.leafGridView();
        const auto& is = gv.indexSet();
        for (const auto& e : elements(gv))
            for (std::size_t k = 0; k < bulkBS; ++k)
                soilRefByCanon[soilCanonicalIndex(e.geometry().center())][k] = solRef[bulkIdx][is.index(e)][k];
    }
    const auto& rootRef = solRef[lowDimIdx]; // indexed by full-network element index == global id

    // ------------------------------------------------------------------
    // 2) parallel: soil split along z; network distributed by the soil partition
    // ------------------------------------------------------------------
    const auto worldComm = typename BulkGrid::Communication(mpiHelper.getCommunicator());
    const auto partitioning = getParam<std::array<int, 3>>(
        "Soil.Grid.Partitioning", std::array<int, 3>{1, 1, worldComm.size()});
    if (partitioning[0] * partitioning[1] * partitioning[2] != worldComm.size())
        DUNE_THROW(Dune::InvalidStateException, "Soil.Grid.Partitioning product must equal the number of ranks");
    Dune::Yasp::FixedSizePartitioning<3> partitioner(partitioning);
    BulkGrid soilParGrid(soilLower, soilUpper, soilCells, std::bitset<3>(0), overlap, worldComm, &partitioner);

    Dune::Communication<Dune::MPIHelper::MPICommunicator> comm(mpiHelper.getCommunicator());
    auto localNetwork = Dune::FoamGridParallel::distributeFromSpatialPartition(
        fullNetwork, soilParGrid.leafGridView(), comm);
    const auto par = localNetwork->parallelData();

    std::vector<double> localRadii(localNetwork->leafGridView().size(0));
    {
        const auto& gv = localNetwork->leafGridView();
        const auto& is = gv.indexSet();
        for (const auto& e : elements(gv))
            localRadii[is.index(e)] = fullRadii[par->elementGlobalId(is.index(e))];
    }

    auto soilParGG = std::make_shared<BulkGridGeometry>(soilParGrid.leafGridView());
    auto rootParGG = std::make_shared<LowDimGridGeometry>(localNetwork->leafGridView());
    auto rootParSpatialParams = std::make_shared<LowDimSpatialParams>(rootParGG, localRadii);
    const auto solPar = solveTransient(soilParGG, rootParGG, rootParSpatialParams);

    // ------------------------------------------------------------------
    // 3) compare owned dofs (all components) and migrated radii against the reference
    // ------------------------------------------------------------------
    std::array<double, bulkBS> bulkDiff{}, bulkRef{};
    std::array<double, lowBS> lowDiff{}, lowRefMax{};
    double maxRadiusDiff = 0.0;
    std::size_t numChecked = 0;

    {
        const auto& gv = soilParGrid.leafGridView();
        const auto& is = gv.indexSet();
        for (const auto& e : elements(gv, Dune::Partitions::interior))
        {
            const auto& ref = soilRefByCanon[soilCanonicalIndex(e.geometry().center())];
            for (std::size_t k = 0; k < bulkBS; ++k)
            {
                bulkDiff[k] = std::max(bulkDiff[k], std::abs(solPar[bulkIdx][is.index(e)][k] - ref[k]));
                bulkRef[k] = std::max(bulkRef[k], std::abs(ref[k]));
            }
            ++numChecked;
        }
    }
    {
        const auto& gv = localNetwork->leafGridView();
        const auto& is = gv.indexSet();
        for (const auto& e : elements(gv, Dune::Partitions::interior))
        {
            const auto gid = par->elementGlobalId(is.index(e));
            for (std::size_t k = 0; k < lowBS; ++k)
            {
                lowDiff[k] = std::max(lowDiff[k], std::abs(solPar[lowDimIdx][is.index(e)][k] - rootRef[gid][k]));
                lowRefMax[k] = std::max(lowRefMax[k], std::abs(rootRef[gid][k]));
            }
            maxRadiusDiff = std::max(maxRadiusDiff, std::abs(localRadii[is.index(e)] - fullRadii[gid]));
            ++numChecked;
        }
    }

    double relErr = 0.0;
    for (std::size_t k = 0; k < bulkBS; ++k)
        relErr = std::max(relErr, comm.max(bulkDiff[k]) / (comm.max(bulkRef[k]) > 0.0 ? comm.max(bulkRef[k]) : 1.0));
    for (std::size_t k = 0; k < lowBS; ++k)
        relErr = std::max(relErr, comm.max(lowDiff[k]) / (comm.max(lowRefMax[k]) > 0.0 ? comm.max(lowRefMax[k]) : 1.0));
    maxRadiusDiff = comm.max(maxRadiusDiff);
    const std::size_t totalChecked = comm.sum(numChecked);
    const bool ok = (relErr < 1e-6) && (maxRadiusDiff == 0.0);

    if (comm.rank() == 0)
    {
        std::cout << "1d3d root-soil 1p2c (transient): compared " << totalChecked << " interior dofs over "
                  << comm.size() << " rank(s); max relative |parallel - sequential| = " << relErr
                  << ", max radius mismatch = " << maxRadiusDiff << "\n";
        std::cout << (ok ? "PASS: parallel transient root-soil solution and migrated radii match sequential\n"
                         : "FAIL: parallel transient root-soil solution/radii differ from sequential\n");
    }

    return ok ? 0 : 1;
}
