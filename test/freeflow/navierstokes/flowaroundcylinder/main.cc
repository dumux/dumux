// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the DFG benchmark of the laminar flow around a cylinder (Navier-Stokes)
 */
#include <config.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/metadata.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_ug.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/mumpssolver.hh>
#include <dumux/linear/stokes_solver.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/assembler.hh>
#include <dumux/multidomain/multistagemultidomainassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

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

    // only the momentum problem has Dirichlet dofs
    if constexpr (requires { momentumProblem->constraints(); })
    {
        for (const auto& constraintData : momentumProblem->constraints())
        {
            const auto& constraintInfo = constraintData.constraintInfo();
            const auto dofIdx = constraintData.dofIndex();
            for (int eqIdx = 0; eqIdx < constraintInfo.size(); ++eqIdx)
                if (constraintInfo.isConstraintEquation(eqIdx))
                    dirichletDofs[momentumIdx][dofIdx][eqIdx] = 1.0;
        }
    }
    else
    {
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
    // run the (in)stationary non-linear problem on this grid
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
    std::cout << "Total number of dofs on rank " << mpiHelper.rank() << ": "
        << massGridGeometry->numDofs() + momentumGridGeometry->numDofs()*MomentumGridGeometry::GridView::dimension
        << std::endl;

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // the instationary cases (2D-2, 2D-3) are only implemented for the CVFE schemes with the new problem interface
    const bool isStationary = !getParam<bool>("Problem.Instationary", false);
#if !NEW_PROBLEM_INTERFACE
    if (!isStationary)
        DUNE_THROW(Dune::NotImplemented, "The instationary cases require the new problem interface (PQ1Bubble/Box)");
#endif

    // compute coupling stencil and afterwards initialize grid variables (need coupling information)
    auto xOld = x;
    if (isStationary)
        couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    else
        couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    // initialize the vtk output module
    constexpr bool isDiamond = MassGridGeometry::discMethod == DiscretizationMethods::fcdiamond;
    const auto mode = isDiamond ? Dune::VTK::nonconforming : Dune::VTK::conforming;
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name(), "", mode);
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    const auto vtkOutputType = isStationary ? Dune::VTK::ascii : Dune::VTK::appendedraw;
    vtkWriter.write(0.0, vtkOutputType);

    // the benchmark indicators are only evaluated with inertia terms (Navier-Stokes)
    const bool evalIndicators = momentumProblem->enableInertiaTerms();
    if (evalIndicators && leafGridView.comm().size() > 1)
        DUNE_THROW(Dune::NotImplemented, "Evaluation of the benchmark indicators in parallel");

    // the benchmark indicators and the grid size are written to a JSON file
    const int dim = MomentumGridGeometry::GridView::dimension;
    MetaData::Collector collector;
    collector["refinement"] = getParam<int>("Grid.Refinement", 0);
    collector["upwindWeight"] = getParam<double>("Flux.UpwindWeight");
    collector["numElements"] = leafGridView.size(0);
    collector["numVertices"] = leafGridView.size(dim);
    collector["numDofs"]["momentum"] = momentumGridGeometry->numDofs();
    collector["numDofs"]["mass"] = massGridGeometry->numDofs();
    collector["numDofs"]["total"] = momentumGridGeometry->numDofs()*dim + massGridGeometry->numDofs();
    const auto indicatorFile = getParam<std::string>("Problem.IndicatorFile", massProblem->name());

    // The iterative Stokes solver is preconditioned for the symmetric Stokes system and is no suitable
    // solver for the Navier-Stokes equations, which are solved with a direct solver
    const bool useStokesSolver = getParam<bool>("LinearSolver.UseIterativeSolver", false);
    if (useStokesSolver && momentumProblem->enableInertiaTerms())
        DUNE_THROW(ParameterException, "LinearSolver.UseIterativeSolver = true (Stokes solver) requires Problem.EnableInertiaTerms = false");

    // create a Newton solver with the chosen linear solver for the given assembler and pass it to apply
    const auto withNewtonSolver = [&](auto assembler, auto&& apply)
    {
        using Assembler = typename decltype(assembler)::element_type;
        const auto applyWithLinearSolver = [&](auto linearSolver)
        {
            using LinearSolver = typename decltype(linearSolver)::element_type;
            using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
            apply(std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager));
        };

        using LATraits = LinearAlgebraTraitsFromAssembler<Assembler>;
        if (useStokesSolver)
        {
            using Matrix = typename Assembler::JacobianMatrix;
            using Vector = typename Assembler::ResidualType;
            auto dDofs = dirichletDofs<Vector>(momentumGridGeometry, massGridGeometry, momentumProblem, momentumIdx, massIdx);
            applyWithLinearSolver(std::make_shared<StokesSolver<Matrix, Vector, MomentumGridGeometry, MassGridGeometry>>(
                momentumGridGeometry, massGridGeometry, dDofs));
        }
#if DUMUX_HAVE_MUMPS && NEW_PROBLEM_INTERFACE
        // MUMPS keeps the symbolic analysis over the many factorizations of the time-dependent cases (same matrix pattern)
        else if (!isStationary)
            applyWithLinearSolver(std::make_shared<DirectSolverMumps<LinearSolverTraits<MomentumGridGeometry>, LATraits>>(
                std::make_tuple(momentumGridGeometry, massGridGeometry)));
#endif
        else
            applyWithLinearSolver(std::make_shared<UMFPackIstlSolver<SeqLinearSolverTraits, LATraits>>());
    };

    if (isStationary)
    {
#if NEW_PROBLEM_INTERFACE
        using Assembler = Experimental::MultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
#else
        using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
#endif
        auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                     std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                     std::make_tuple(momentumGridVariables, massGridVariables),
                                                     couplingManager);
        withNewtonSolver(assembler, [&](auto nonLinearSolver) { nonLinearSolver->solve(x); });
        vtkWriter.write(1.0);
        couplingManager->updateSolution(x);

        // evaluate benchmark indicators (the reference values are only valid for case 2D-1 at Re = 20)
        if (evalIndicators)
        {
            static constexpr double cDragReference = 5.57953523384;
            static constexpr double cLiftReference = 0.010618948146;
            static constexpr double pDiffReference = 0.11752016697;
            const auto [cDrag, cLift] = momentumProblem->evalDragAndLiftCoefficient(*momentumGridVariables, x[momentumIdx]);
            const auto pDiff = massProblem->evalPressureDifference(*massGridVariables, x[massIdx]);
            std::cout << "cDrag: " << cDrag << " (reference: " << cDragReference << ")\n"
                      << "cLift: " << cLift << " (reference: " << cLiftReference << ")\n"
                      << "pDiff: " << pDiff << " (reference: " << pDiffReference << ")" << std::endl;

            collector["cDrag"] = cDrag;
            collector["cLift"] = cLift;
            collector["pDiff"] = pDiff;
            MetaData::writeJsonFile(collector, indicatorFile);

            // To be close to the reference values the grid needs to be refined
            if (getParam<bool>("Problem.CheckIndicators", false))
            {
                using std::abs;
                if (abs(cDrag - cDragReference) > 0.002)
                    DUNE_THROW(Dune::Exception,
                        "Deviation from reference too large for drag coefficient: "
                        << cDrag << " (ref: " << cDragReference << ")"
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
    }
#if NEW_PROBLEM_INTERFACE
    else
    {
        // DIRK schemes do not weight the pressure of the old time level, so the indicators can be evaluated at
        // the new time level; with Crank-Nicolson the pressure oscillates (see README, evaluation at the time levels)
        using Scalar = typename Traits::Scalar;
        const auto scheme = getParam<std::string>("TimeLoop.Scheme", "DIRK2");
        const auto method = [&]() -> std::shared_ptr<Experimental::MultiStageMethod<Scalar>>
        {
            if (scheme == "DIRK2")
                return std::make_shared<Experimental::MultiStage::DIRKSecondOrderAlexander<Scalar>>();
            else if (scheme == "DIRK3")
                return std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<Scalar>>();
            else if (scheme == "ImplicitEuler")
                return std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();
            else if (scheme == "CrankNicolson")
                return std::make_shared<Experimental::MultiStage::Theta<Scalar>>(0.5);
            else
                DUNE_THROW(ParameterException, "Unknown time stepping scheme " << scheme
                           << ". Possible values are DIRK2, DIRK3, ImplicitEuler and CrankNicolson.");
        }();
        collector["timeSteppingScheme"] = method->name();
        if (scheme == "CrankNicolson" && evalIndicators)
            std::cout << "Warning: with Crank-Nicolson, the pressure at the time levels oscillates and the"
                      << " pressure-dependent benchmark indicators are not meaningful" << std::endl;
        using Assembler = Experimental::MultiStageMultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
        auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                     std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                     std::make_tuple(momentumGridVariables, massGridVariables),
                                                     couplingManager, method, xOld);

        // Fixed time step sizes: Dt until TEnd (case 2D-3), or a coarse step Dt1 until T1 followed
        // by a fine step Dt2 until T2 (case 2D-2, the coarse phase develops the periodic flow)
        const bool twoPhases = hasParam("TimeLoop.T1");
        const auto tEnd = getParam<Scalar>(twoPhases ? "TimeLoop.T2" : "TimeLoop.TEnd");
        const auto dt = getParam<Scalar>(twoPhases ? "TimeLoop.Dt2" : "TimeLoop.Dt");
        const auto startupTEnd = twoPhases ? getParam<Scalar>("TimeLoop.T1") : 0.0;
        const auto startupDt = twoPhases ? getParam<Scalar>("TimeLoop.Dt1") : dt;
        const auto fixedTimeStepSize = [&](const Scalar time)
        {
            // the last coarse step ends exactly at T1 (unless it would be shorter than half a step)
            using std::min;
            return time < startupTEnd - 0.5*startupDt ? min(startupDt, startupTEnd - time) : dt;
        };

        auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, fixedTimeStepSize(0.0), tEnd);
        using std::max;
        timeLoop->setMaxTimeStepSize(max(dt, startupDt));
        // VTK output at fixed times (every step if not set), independent of the time step size
        const auto vtkOutputTimeInterval = getParam<Scalar>("Vtk.OutputTimeInterval", 0.0);
        Scalar nextVtkOutputTime = vtkOutputTimeInterval;
        const auto isVtkOutputTime = [&](const Scalar time)
        {
            if (vtkOutputTimeInterval <= 0.0)
                return true;
            const Scalar eps = 1e-6*timeLoop->timeStepSize();
            if (time < nextVtkOutputTime - eps)
                return false;
            while (nextVtkOutputTime <= time + eps)
                nextVtkOutputTime += vtkOutputTimeInterval;
            return true;
        };

        // the indicators are appended to a CSV file at every time step (flushed), so they survive an aborted
        // run (e.g. killed by the wall-time limit of a cluster); the JSON file is only written at the end
        std::ofstream indicatorCsv;
        if (evalIndicators)
        {
            indicatorCsv.open(indicatorFile + ".csv");
            indicatorCsv << "time,cDrag,cLift,pDiff" << std::endl;
            indicatorCsv << std::setprecision(std::numeric_limits<Scalar>::max_digits10);
        }

        const auto collectIndicators = [&](const Scalar time)
        {
            // the traction on the cylinder uses the pressure of the coupling manager
            couplingManager->updateSolution(x);
            const auto [cDrag, cLift] = momentumProblem->evalDragAndLiftCoefficient(*momentumGridVariables, x[momentumIdx], false);
            const auto pDiff = massProblem->evalPressureDifference(*massGridVariables, x[massIdx]);

            std::cout << "Benchmark indicators at t = " << time << ": cDrag = " << cDrag
                      << ", cLift = " << cLift << ", pDiff = " << pDiff << std::endl;
            collector["time"].push_back(time);
            collector["cDrag"].push_back(cDrag);
            collector["cLift"].push_back(cLift);
            collector["pDiff"].push_back(pDiff);
            indicatorCsv << time << "," << cDrag << "," << cLift << "," << pDiff << std::endl;
        };

        // the indicators are written at the end and if the simulation fails
        const auto writeIndicators = [&](const bool finished, const std::string& message)
        {
            if (!evalIndicators)
                return;
            collector["finished"] = finished;
            if (!finished)
                collector["error"] = message;
            MetaData::writeJsonFile(collector, indicatorFile);
        };

        try
        {
            withNewtonSolver(assembler, [&](auto nonLinearSolver)
            {
                using NewtonSolver = typename decltype(nonLinearSolver)::element_type;
                Experimental::MultiStageTimeStepper<NewtonSolver> timeStepper(nonLinearSolver, method);

                timeLoop->start(); do
                {
                    // solve the time step (with a reduced time step size if Newton's method fails)
                    timeStepper.step(x, *timeLoop);

                    if (evalIndicators)
                        collectIndicators(timeLoop->time() + timeLoop->timeStepSize());

                    xOld = x;
                    momentumGridVariables->advanceTimeStep();
                    massGridVariables->advanceTimeStep();
                    timeLoop->advanceTimeStep();

                    if (isVtkOutputTime(timeLoop->time()) || timeLoop->finished())
                        vtkWriter.write(timeLoop->time(), vtkOutputType);

                    timeLoop->reportTimeStep();

                    // continue with the fixed time step size (also after a reduction on failure)
                    timeLoop->setTimeStepSize(fixedTimeStepSize(timeLoop->time()));

                } while (!timeLoop->finished());

                timeLoop->finalize(leafGridView.comm());
            });
        }
        catch (const Dune::Exception& e)
        {
            writeIndicators(false, e.what());
            throw;
        }
        catch (const std::exception& e)
        {
            writeIndicators(false, e.what());
            throw;
        }

        writeIndicators(true, "");
    }
#endif

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
