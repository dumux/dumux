// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief Production driver for the MPI-parallel transient root-soil benchmark (1p2c Richards /
 *        1pnc root) on a realistic DGF network, overlapping (Strategy A) decomposition with the
 *        parallel multidomain MUMPS backend.
 *
 * Unlike the regression test (which also solves a replicated reference and compares), this runs
 * only the distributed problem over the real time span with episode-based VTK output, a parallel
 * global mass balance, and per-phase wall-clock timing for scaling studies.
 */
#include <config.h>

#include <array>
#include <bitset>
#include <iostream>
#include <iomanip>
#include <memory>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
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
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/mumpssolver.hh>
#if defined(USE_RAS_SOLVER)
// Block-diagonal (Restricted Additive Schwarz) parallel solver: exact local MUMPS block solve
// per rank inside a parallel BiCGSTAB, an alternative to the global parallel direct factorization
// for the scaling study.
#include <dumux/linear/mumpsrassolver.hh>
#elif defined(USE_BLOCKDIAGILU_SOLVER)
// The sequential field-block-diagonal ILU0-preconditioned BiCGSTAB used by the original
// (non-parallel) embedded 1p2c/richards2c test. Sequential only (no MPI); for baseline comparison.
#include <dumux/linear/seqsolverbackend.hh>
#endif

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include "properties.hh"

namespace Dumux {

//! global tracer mass over owned (interior) cells, reduced over all ranks
template<class Problem, class SolutionVector, class GridVariables, class Comm>
double globalTracerMass(const Problem& problem, const SolutionVector& sol, const GridVariables& gridVars, const Comm& comm)
{
    static constexpr int liquidPhaseIdx = Problem::Indices::liquidPhaseIdx;
    static constexpr int transportCompIdx = Problem::Indices::transportCompIdx;
    double mass = 0.0;
    const auto& gg = problem.gridGeometry();
    auto fvGeometry = localView(gg);
    auto elemVolVars = localView(gridVars.curGridVolVars());
    for (const auto& element : elements(gg.gridView(), Dune::Partitions::interior))
    {
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, sol);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& vv = elemVolVars[scv];
            mass += vv.massFraction(liquidPhaseIdx, transportCompIdx)*vv.density(liquidPhaseIdx)
                    *scv.volume()*vv.porosity()*vv.saturation(liquidPhaseIdx)*vv.extrusionFactor();
        }
    }
    return comm.sum(mass);
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    Parameters::init(argc, argv);

    using BulkTypeTag = Properties::TTag::Soil;
    using LowDimTypeTag = Properties::TTag::Root;

    using BulkGrid = GetPropType<BulkTypeTag, Properties::Grid>;
    using LowDimGrid = GetPropType<LowDimTypeTag, Properties::Grid>;
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
#elif defined(USE_BLOCKDIAGILU_SOLVER)
    using DirectSolver = BlockDiagILU0BiCGSTABSolver;
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
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize", tEnd);
    const auto episodeLength = getParam<Scalar>("TimeLoop.EpisodeLength", tEnd);
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    const bool outputVtk = getParam<bool>("Problem.EnableVtkOutput", false);

    Dune::Communication<Dune::MPIHelper::MPICommunicator> comm(mpiHelper.getCommunicator());
    const int rank = comm.rank();

    // ---- build the distributed grids (soil split; network nested in the soil partition) ----
    const auto worldComm = typename BulkGrid::Communication(mpiHelper.getCommunicator());
    const auto partitioning = getParam<std::array<int, 3>>(
        "Soil.Grid.Partitioning", std::array<int, 3>{1, 1, worldComm.size()});
    if (partitioning[0] * partitioning[1] * partitioning[2] != worldComm.size())
        DUNE_THROW(Dune::InvalidStateException, "Soil.Grid.Partitioning product must equal the number of ranks");
    Dune::Yasp::FixedSizePartitioning<3> partitioner(partitioning);
    BulkGrid soilGrid(soilLower, soilUpper, soilCells, std::bitset<3>(0), overlap, worldComm, &partitioner);

    Dune::GridPtr<LowDimGrid> networkPtr(networkFile, mpiHelper.getLocalCommunicator());
    auto& fullNetwork = *networkPtr;
    auto localNetwork = Dune::FoamGridParallel::distributeFromSpatialPartition(
        fullNetwork, soilGrid.leafGridView(), comm);
    const auto par = localNetwork->parallelData();

    // migrate per-element radius (from the DGF surface parameter) by global id
    std::vector<double> fullRadii(fullNetwork.leafGridView().size(0));
    for (const auto& e : elements(fullNetwork.leafGridView()))
        fullRadii[fullNetwork.leafGridView().indexSet().index(e)]
            = networkPtr.parameters(e)[surfaceIdx] / e.geometry().volume() / (2.0 * M_PI);
    std::vector<double> localRadii(localNetwork->leafGridView().size(0));
    for (const auto& e : elements(localNetwork->leafGridView()))
        localRadii[localNetwork->leafGridView().indexSet().index(e)]
            = fullRadii[par->elementGlobalId(localNetwork->leafGridView().indexSet().index(e))];

    auto soilGG = std::make_shared<BulkGridGeometry>(soilGrid.leafGridView());
    auto rootGG = std::make_shared<LowDimGridGeometry>(localNetwork->leafGridView());

    if (rank == 0)
        std::cout << "Distributed problem: " << comm.size() << " rank(s), "
                  << "soil " << soilCells[0] << "x" << soilCells[1] << "x" << soilCells[2]
                  << ", network " << fullNetwork.leafGridView().size(0) << " elements\n";

    // ---- set up the coupled problem ----
    auto couplingManager = std::make_shared<CouplingManager>(soilGG, rootGG);
    auto soilProblem = std::make_shared<BulkProblem>(soilGG, couplingManager);
    auto rootSpatialParams = std::make_shared<LowDimSpatialParams>(rootGG, localRadii);
    auto rootProblem = std::make_shared<LowDimProblem>(rootGG, rootSpatialParams, couplingManager);

    SolutionVector sol;
    sol[bulkIdx].resize(soilGG->numDofs());
    sol[lowDimIdx].resize(rootGG->numDofs());
    soilProblem->applyInitialSolution(sol[bulkIdx]);
    rootProblem->applyInitialSolution(sol[lowDimIdx]);
    auto oldSol = sol;

    couplingManager->init(soilProblem, rootProblem, sol);
    soilProblem->computePointSourceMap();
    rootProblem->computePointSourceMap();

    auto soilGridVariables = std::make_shared<BulkGridVariables>(soilProblem, soilGG);
    soilGridVariables->init(sol[bulkIdx]);
    auto rootGridVariables = std::make_shared<LowDimGridVariables>(rootProblem, rootGG);
    rootGridVariables->init(sol[lowDimIdx]);

    // ---- vtk output (parallel: each rank writes its .vtu, rank 0 writes the .pvtu) ----
    using BulkSolutionVector = std::decay_t<decltype(sol[bulkIdx])>;
    using LowDimSolutionVector = std::decay_t<decltype(sol[lowDimIdx])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> soilVtk(*soilGridVariables, sol[bulkIdx], soilProblem->name());
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> rootVtk(*rootGridVariables, sol[lowDimIdx], rootProblem->name());
    if (outputVtk)
    {
        GetPropType<BulkTypeTag, Properties::IOFields>::initOutputModule(soilVtk);
        GetPropType<LowDimTypeTag, Properties::IOFields>::initOutputModule(rootVtk);
        rootProblem->addVtkOutputFields(rootVtk);
        soilVtk.write(0.0);
        rootVtk.write(0.0);
    }

    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    timeLoop->setPeriodicCheckPoint(episodeLength);

    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(soilProblem, rootProblem),
        std::make_tuple(soilGG, rootGG),
        std::make_tuple(soilGridVariables, rootGridVariables),
        couplingManager, timeLoop, oldSol);

#if defined(USE_BLOCKDIAGILU_SOLVER)
    if (comm.size() > 1)
        DUNE_THROW(Dune::NotImplemented, "BlockDiagILU0BiCGSTABSolver is sequential only; run with one rank");
    auto linearSolver = std::make_shared<DirectSolver>();
#else
    std::shared_ptr<DirectSolver> linearSolver = (comm.size() > 1)
        ? std::make_shared<DirectSolver>(std::make_tuple(soilGG, rootGG))
        : std::make_shared<DirectSolver>();
#endif
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager, comm);

    const double initialMass = globalTracerMass(*soilProblem, sol[bulkIdx], *soilGridVariables, comm)
                             + globalTracerMass(*rootProblem, sol[lowDimIdx], *rootGridVariables, comm);
    if (rank == 0)
        std::cout << "Initial tracer mass: " << initialMass*1e12 << " ng\n";

    // ---- time loop (timed for scaling) ----
    Dune::Timer wallClock;
    timeLoop->start();
    while (!timeLoop->finished())
    {
        nonLinearSolver.solve(sol, *timeLoop);
        oldSol = sol;
        soilGridVariables->advanceTimeStep();
        rootGridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        if (outputVtk && (timeLoop->isCheckPoint() || timeLoop->finished()))
        {
            soilVtk.write(timeLoop->time());
            rootVtk.write(timeLoop->time());
        }

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    }
    const double loopTime = wallClock.elapsed();

    const double finalMass = globalTracerMass(*soilProblem, sol[bulkIdx], *soilGridVariables, comm)
                           + globalTracerMass(*rootProblem, sol[lowDimIdx], *rootGridVariables, comm);

    if (rank == 0)
    {
        std::cout << "Final tracer mass: " << finalMass*1e12 << " ng\n";
        std::cout << "[timing] ranks=" << comm.size()
                  << " soilCells=" << std::size_t(soilCells[0])*soilCells[1]*soilCells[2]
                  << " networkElems=" << fullNetwork.leafGridView().size(0)
                  << " timeLoopWallTime=" << std::fixed << std::setprecision(3) << loopTime << " s"
                  << " steps=" << timeLoop->timeStepIndex() << "\n";
    }

    return 0;
}
