// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/common/version.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/stokes_solver.hh>
#include <dumux/linear/navierstokes_simple_solver.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
// The Amesos2 (trilinossolvers.hh) and direct MUMPS (mumpssolver.hh) backends are mutually
// exclusive in one translation unit (Amesos2 drags dmumps_c.h into its own namespace, which
// shadows the global MUMPS struct). Select one at compile time. Both expose the same direct
// solver interface (two-grid-geometry constructor in parallel, no-arg constructor sequentially).
#if defined(USE_MUMPS_SOLVER)
#include <dumux/linear/mumpssolver.hh>
#define DUMUX_TEST_DIRECT_SOLVER ::Dumux::DirectSolverMumps
#define DUMUX_TEST_HAVE_DIRECT DUMUX_HAVE_MUMPS
#else
#if DUMUX_HAVE_TRILINOS
#include <dumux/linear/trilinossolvers.hh>
#endif
#define DUMUX_TEST_DIRECT_SOLVER ::Dumux::DirectSolverAmesos2
#define DUMUX_TEST_HAVE_DIRECT DUMUX_HAVE_TRILINOS
#endif

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/momentum/velocityoutput.hh>

#include "properties.hh"

template<class Vector, class MomGG, class MassGG, class MomP, class MomIdx, class MassIdx>
auto dirichletDofs(std::shared_ptr<MomGG> momentumGridGeometry,
                   std::shared_ptr<MassGG> massGridGeometry,
                   std::shared_ptr<MomP> momentumProblem,
                   MomIdx momentumIdx, MassIdx massIdx)
{
    Vector dirichletDofs;
    dirichletDofs[momentumIdx].resize(momentumGridGeometry->numDofs());
    dirichletDofs[massIdx].resize(massGridGeometry->numDofs());
    dirichletDofs = 0.0;

    auto fvGeometry = localView(*momentumGridGeometry);
    for (const auto& element : elements(momentumGridGeometry->gridView()))
    {
        fvGeometry.bind(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            if (momentumGridGeometry->dofOnBoundary(scv.dofIndex()))
            {
                const auto bcTypes = momentumProblem->boundaryTypes(element, scv);
                for (int i = 0; i < bcTypes.size(); ++i)
                    if (bcTypes.isDirichlet(i))
                        dirichletDofs[momentumIdx][scv.dofIndex()][i] = 1.0;
            }
        }
    }

    return dirichletDofs;
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::TYPETAG_MOMENTUM;
    using MassTypeTag = Properties::TTag::TYPETAG_MASS;

    // initialize MPI, finalize is done automatically on exit
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    Dune::Timer timer;

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // create a grid
    using Grid = GetPropType<MassTypeTag, Properties::Grid>;
    Dumux::GridManager<Grid> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    std::cout << "Grid overlap size: " << leafGridView.overlapSize(0) << std::endl;
    std::cout << "Grid ghost size: " << leafGridView.ghostSize(0) << std::endl;

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // the coupling manager
    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());
    x = 0.0;
    auto xOld = x;
    std::cout << "Total number of dofs on rank " << mpiHelper.rank() << ": "
        << massGridGeometry->numDofs() + momentumGridGeometry->numDofs()*MomentumGridGeometry::GridView::dimension
        << std::endl;

    // optional instationary (time-dependent) Navier-Stokes: build a time loop when requested
    using Scalar = typename Traits::Scalar;
    const bool instationary = getParam<bool>("Problem.Instationary", false);
    std::shared_ptr<TimeLoop<Scalar>> timeLoop;
    if (instationary)
    {
        const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
        const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
        const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize", tEnd);
        timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
        timeLoop->setMaxTimeStepSize(maxDt);
    }

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // compute coupling stencil and afterwards initialize grid variables (need coupling information)
    if (instationary)
        couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);
    else
        couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    std::shared_ptr<Assembler> assembler;
    if (instationary)
        assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                std::make_tuple(momentumGridVariables, massGridVariables),
                                                couplingManager, timeLoop, xOld);
    else
        assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                std::make_tuple(momentumGridVariables, massGridVariables),
                                                couplingManager);

    // initialize the vtk output module
    constexpr bool isDiamond = MassGridGeometry::discMethod == DiscretizationMethods::fcdiamond;
    const auto mode = isDiamond ? Dune::VTK::nonconforming : Dune::VTK::conforming;
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name(), "", mode);
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.write(0.0);

    // Benchmark-indicator time series. The DFG 2D-2 (Re=100) case develops an oscillatory
    // steady state (Karman vortex street), so the drag/lift coefficients and the pressure
    // difference become time-periodic -- see benchmarkreferences.md (dfg_benchmark2_re100).
    // We log them at every time step (rank 0) so the onset of oscillation is observable.
    const bool writeIndicators = instationary && momentumProblem->enableInertiaTerms();
    std::ofstream benchmarkCsv;
    if (writeIndicators && mpiHelper.rank() == 0)
    {
        benchmarkCsv.open(massProblem->name() + "_benchmark_indicators.csv");
        benchmarkCsv << "time,cDrag,cLift,pDiff\n" << std::setprecision(10);
    }

    // run either a single stationary solve or a time loop (instationary)
    auto runSimulation = [&](auto& nonLinearSolver)
    {
        if (instationary)
        {
            timeLoop->start(); do
            {
                // solve the non-linear system with time step control
                nonLinearSolver.solve(x, *timeLoop);

                // make the new solution the old solution
                xOld = x;
                momentumGridVariables->advanceTimeStep();
                massGridVariables->advanceTimeStep();

                // advance the time loop to the next step and write output
                timeLoop->advanceTimeStep();
                vtkWriter.write(timeLoop->time());
                timeLoop->reportTimeStep();
                timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

                // benchmark indicators for this time step (drag/lift/pressure difference).
                // eval* do collective MPI reductions, so all ranks must call them.
                if (writeIndicators)
                {
                    couplingManager->updateSolution(x);
                    const auto [cDrag, cLift] = momentumProblem->evalDragAndLiftCoefficient(*momentumGridVariables, x[momentumIdx]);
                    const auto pDiff = massProblem->evalPressureDifference(*massGridVariables, x[massIdx]);
                    if (mpiHelper.rank() == 0)
                    {
                        std::cout << "Benchmark indicators @ t = " << timeLoop->time()
                                  << ":  cDrag = " << cDrag << ",  cLift = " << cLift
                                  << ",  pDiff = " << pDiff << std::endl;
                        benchmarkCsv << timeLoop->time() << "," << cDrag << ","
                                     << cLift << "," << pDiff << "\n";
                        benchmarkCsv.flush();
                    }
                }

            } while (!timeLoop->finished());

            timeLoop->finalize(leafGridView.comm());
        }
        else
            nonLinearSolver.solve(x);
    };

    // the linearize and solve
    if (getParam<bool>("LinearSolver.UseSimpleSolver", false))
    {
        using Matrix = typename Assembler::JacobianMatrix;
        using Vector = typename Assembler::ResidualType;
        using LinearSolver = NavierStokesSimpleSolver<Matrix, Vector, MomentumGridGeometry, MassGridGeometry>;
        auto dDofs = dirichletDofs<Vector>(momentumGridGeometry, massGridGeometry, momentumProblem, momentumIdx, massIdx);
        auto linearSolver = std::make_shared<LinearSolver>(momentumGridGeometry, massGridGeometry, dDofs);
        using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
        NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
        runSimulation(nonLinearSolver);
    }
    else if (getParam<bool>("LinearSolver.UseIterativeSolver", false))
    {
        using Matrix = typename Assembler::JacobianMatrix;
        using Vector = typename Assembler::ResidualType;
        using LinearSolver = StokesSolver<Matrix, Vector, MomentumGridGeometry, MassGridGeometry>;
        auto dDofs = dirichletDofs<Vector>(momentumGridGeometry, massGridGeometry, momentumProblem, momentumIdx, massIdx);
        auto linearSolver = std::make_shared<LinearSolver>(momentumGridGeometry, massGridGeometry, dDofs);
        using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
        NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
        runSimulation(nonLinearSolver);
    }
    else if (getParam<bool>("LinearSolver.UseTrilinos", false) || getParam<bool>("LinearSolver.UseMumps", false))
    {
#if DUMUX_TEST_HAVE_DIRECT
        using LATraits = LinearAlgebraTraitsFromAssembler<Assembler>;
        using DirectSolver = DUMUX_TEST_DIRECT_SOLVER<SeqLinearSolverTraits, LATraits>;
        if (mpiHelper.size() > 1)
        {
            // Parallel: pass the subdomain grid geometries as a tuple to build combined DOF maps
            auto linearSolver = std::make_shared<DirectSolver>(std::make_tuple(momentumGridGeometry, massGridGeometry));
            using NewtonSolver = MultiDomainNewtonSolver<Assembler, DirectSolver, CouplingManager>;
            NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
            runSimulation(nonLinearSolver);
        }
        else
        {
            // Sequential: use the no-arg constructor with lazy DOF map init
            auto linearSolver = std::make_shared<DirectSolver>();
            using NewtonSolver = MultiDomainNewtonSolver<Assembler, DirectSolver, CouplingManager>;
            NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
            runSimulation(nonLinearSolver);
        }
#else
        DUNE_THROW(Dune::NotImplemented, "Direct solver backend not available");
#endif
    }
    else
    {
        using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
        auto linearSolver = std::make_shared<LinearSolver>();
        using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
        NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);
        runSimulation(nonLinearSolver);
    }

    // write vtk output (instationary already writes inside the time loop)
    if (!instationary)
        vtkWriter.write(1.0);

    // update coupling manager for output
    couplingManager->updateSolution(x);

    // evaluate benchmark indicators (drag/lift/pressure-difference; computed in parallel too).
    // For the instationary case the indicators are logged per time step in the time loop above.
    if (!instationary && momentumProblem->enableInertiaTerms())
    {
        static constexpr double cDrafReference = 5.57953523384;
        static constexpr double cLiftReference = 0.010618948146;
        static constexpr double pDiffReference = 0.11752016697;
        const auto [cDrag, cLift] = momentumProblem->evalDragAndLiftCoefficient(*momentumGridVariables, x[momentumIdx]);
        std::cout << "cDrag: " << cDrag
                << " (reference: " << cDrafReference << ")"
                << "\n"
                << "cLift: " << cLift
                << " (reference: " << cLiftReference << ")"
                << std::endl;

        const auto pDiff = massProblem->evalPressureDifference(*massGridVariables, x[massIdx]);
        std::cout << "pDiff: " << pDiff
                << " (reference: " << pDiffReference << ")"
                << std::endl;

        // To be close to the reference values the grid needs to be refined
        if (getParam<bool>("Problem.CheckIndicators", false))
        {
            using std::abs;
            if (abs(cDrag - cDrafReference) > 0.002)
                DUNE_THROW(Dune::Exception,
                    "Deviation from reference too large for drag coefficient: "
                    << cDrag << " (ref: " << cDrafReference << ")"
                );
            if (abs(cLift - cLiftReference) > 0.002)
                DUNE_THROW(Dune::Exception,
                    "Deviation from reference too large for lift coefficient: "
                    << cLift << " (ref: " << cLiftReference << ")"
                );
            if (abs(pDiff - pDiffReference) > 0.003)
                DUNE_THROW(Dune::Exception,
                    "Deviation from reference too large for pressure difference: "
                    << pDiff << " (ref: " << pDiffReference << ")"
                );
        }
    }

    timer.stop();
    const auto& comm = leafGridView.comm();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
