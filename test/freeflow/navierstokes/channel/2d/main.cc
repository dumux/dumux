// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model.
 */

#include <config.h>

#ifndef NONISOTHERMAL
#define NONISOTHERMAL 0
#endif

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/stokes_solver.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/freeflow/navierstokes/momentum/velocityoutput.hh>

#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>
#include <test/freeflow/navierstokes/errors.hh>

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
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary() || !scvf.isFrontal())
                continue;
            const auto& scv = fvGeometry.scv(scvf.insideScvIdx());
            if (scv.boundary())
            {
                const auto bcTypes = momentumProblem->boundaryTypes(element, scvf);
                if (bcTypes.isDirichlet(scv.dofAxis()))
                    dirichletDofs[momentumIdx][scv.dofIndex()][0] = 1.0;
            }
        }
    }

    return dirichletDofs;
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::ChannelTestMomentum;
    using MassTypeTag = Properties::TTag::ChannelTestMass;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<MassTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);

    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // the coupling manager
    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);

    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // get some time loop parameters
    using Scalar = GetPropType<MassTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());
    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);
    auto xOld = x;

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    massProblem->setTime(timeLoop->time());
    momentumProblem->setTime(timeLoop->time());

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);

    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // initialize the coupling stencils
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);
    // initializing the gridvariables requires the coupling manager to be set up
    momentumGridVariables->init(x[momentumIdx]);
    massGridVariables->init(x[massIdx]);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());

    using NSAnalyticSol = NavierStokesTest::AnalyticalSolutionVectors<MomentumProblem, MassProblem>;
    std::unique_ptr<NSAnalyticSol> analyticalSolVectors;
    if (momentumProblem->hasAnalyticalSolution())
    {
        analyticalSolVectors = std::make_unique<NSAnalyticSol>(momentumProblem, massProblem);
        vtkWriter.addField(analyticalSolVectors->analyticalPressureSolution(), "pressureExact");
        vtkWriter.addField(analyticalSolVectors->analyticalVelocitySolution(), "velocityExact");
        // vtkWriter.addFaceField(analyticalSolVectors.analyticalVelocitySolutionOnFace(), "faceVelocityExact");
    }

    vtkWriter.write(restartTime);

    // the assembler with time loop for instationary problem
    const bool isStationary = getParam<bool>("Problem.IsStationary", false);
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = isStationary ? std::make_shared<Assembler>(
            std::make_tuple(momentumProblem, massProblem),
            std::make_tuple(momentumGridGeometry, massGridGeometry),
            std::make_tuple(momentumGridVariables, massGridVariables),
            couplingManager
        ) : std::make_shared<Assembler>(
            std::make_tuple(momentumProblem, massProblem),
            std::make_tuple(momentumGridGeometry, massGridGeometry),
            std::make_tuple(momentumGridVariables, massGridVariables),
            couplingManager,
            timeLoop, xOld
        );

#if !NONISOTHERMAL
    const bool useIterativeSolver = getParam<bool>("LinearSolver.UseIterativeSolver", false);
    using Matrix = typename Assembler::JacobianMatrix;
    using Vector = typename Assembler::ResidualType;
    using IterativeLinearSolver = StokesSolver<Matrix, Vector, MomentumGridGeometry, MassGridGeometry>;
    auto dDofs = dirichletDofs<Vector>(momentumGridGeometry, massGridGeometry, momentumProblem, momentumIdx, massIdx);
    auto iterativeLinearSolver = std::make_shared<IterativeLinearSolver>(momentumGridGeometry, massGridGeometry, dDofs);
    using NewtonSolverWithIterativeSolver = MultiDomainNewtonSolver<Assembler, IterativeLinearSolver, CouplingManager>;
    NewtonSolverWithIterativeSolver nonLinearSolverWithIterativeSolver(assembler, iterativeLinearSolver, couplingManager);
#endif

    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    if (isStationary)
    {
        // solve the non-linear system with time step control
#if !NONISOTHERMAL
        if (useIterativeSolver)
            nonLinearSolverWithIterativeSolver.solve(x);
        else
#endif
            nonLinearSolver.solve(x);

        // print discrete L2 and Linfity errors
        if (momentumProblem->hasAnalyticalSolution())
        {
            // print discrete L2 and Linfity errors
            const bool printErrors = getParam<bool>("Problem.PrintErrors", false);
            const bool printConvergenceTestFile = getParam<bool>("Problem.PrintConvergenceTestFile", false);

            if (printErrors || printConvergenceTestFile)
            {
                NavierStokesTest::Errors errors(momentumProblem, massProblem, x);
                NavierStokesTest::ErrorCSVWriter(
                    momentumProblem, massProblem, std::to_string(x[massIdx].size())
                ).printErrors(errors);

                if (printConvergenceTestFile)
                    NavierStokesTest::convergenceTestAppendErrors(momentumProblem, massProblem, errors);
            }

        }
        // write vtk output
        vtkWriter.write(1.0);
    }
    else
    {
        if (getParam<Scalar>("Problem.InletVelocity") > 1e-6)
            timeLoop->setCheckPoint({200.0, 210.0});

        // time loop
        timeLoop->start(); do
        {
            // solve the non-linear system with time step control
#if !NONISOTHERMAL
            if (useIterativeSolver)
                nonLinearSolverWithIterativeSolver.solve(x, *timeLoop);
            else
#endif
                nonLinearSolver.solve(x, *timeLoop);

            // make the new solution the old solution
            xOld = x;
            momentumGridVariables->advanceTimeStep();
            massGridVariables->advanceTimeStep();

            // advance the time loop to the next step
            timeLoop->advanceTimeStep();

            // write vtk output
            vtkWriter.write(timeLoop->time());

            // report statistics of this time step
            timeLoop->reportTimeStep();

            // set new dt as suggested by newton solver
#if !NONISOTHERMAL
            if (useIterativeSolver)
                timeLoop->setTimeStepSize(nonLinearSolverWithIterativeSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
            else
#endif
                timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

            // update time
            massProblem->setTime(timeLoop->time());
            momentumProblem->setTime(timeLoop->time());

        } while (!timeLoop->finished());

        timeLoop->finalize(leafGridView.comm());
    }

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
