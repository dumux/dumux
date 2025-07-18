// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoundaryTests
 * \brief A test problem for the coupled Stokes/Darcy problem (1p).
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/discretization/method.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/freeflow/navierstokes/momentum/velocityoutput.hh>

#include "../properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using FreeFlowMomentumTypeTag = Properties::TTag::FreeFlowOnePMomentum;
    using FreeFlowMassTypeTag = Properties::TTag::FreeFlowOnePMass;
    using DarcyTypeTag = Properties::TTag::DarcyOneP;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using DarcyGridManager = Dumux::GridManager<GetPropType<DarcyTypeTag, Properties::Grid>>;
    DarcyGridManager darcyGridManager;
    darcyGridManager.init("Darcy"); // pass parameter group

    using FreeFlowGridManager = Dumux::GridManager<GetPropType<FreeFlowMomentumTypeTag, Properties::Grid>>;
    FreeFlowGridManager freeFlowGridManager;
    freeFlowGridManager.init("FreeFlow"); // pass parameter group

    // we compute on the leaf grid view
    const auto& darcyGridView = darcyGridManager.grid().leafGridView();
    const auto& freeFlowGridView = freeFlowGridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FreeFlowMomentumGridGeometry = GetPropType<FreeFlowMomentumTypeTag, Properties::GridGeometry>;
    auto freeFlowMomentumGridGeometry = std::make_shared<FreeFlowMomentumGridGeometry>(freeFlowGridView);
    using FreeFlowMassGridGeometry = GetPropType<FreeFlowMassTypeTag, Properties::GridGeometry>;
    auto freeFlowMassGridGeometry = std::make_shared<FreeFlowMassGridGeometry>(freeFlowGridView);
    using DarcyGridGeometry = GetPropType<DarcyTypeTag, Properties::GridGeometry>;
    auto darcyGridGeometry = std::make_shared<DarcyGridGeometry>(darcyGridView);

    using Traits = MultiDomainTraits<FreeFlowMomentumTypeTag, FreeFlowMassTypeTag, DarcyTypeTag>;
    using CouplingManager = FreeFlowPorousMediumCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the indices
    constexpr auto freeFlowMomentumIndex = CouplingManager::freeFlowMomentumIndex;
    constexpr auto freeFlowMassIndex = CouplingManager::freeFlowMassIndex;
    constexpr auto porousMediumIndex = CouplingManager::porousMediumIndex;

    // the problem (initial and boundary conditions)
    using FreeFlowMomentumProblem = GetPropType<FreeFlowMomentumTypeTag, Properties::Problem>;
    auto freeFlowMomentumProblem = std::make_shared<FreeFlowMomentumProblem>(freeFlowMomentumGridGeometry, couplingManager);
    using FreeFlowMassProblem = GetPropType<FreeFlowMassTypeTag, Properties::Problem>;
    auto freeFlowMassProblem = std::make_shared<FreeFlowMassProblem>(freeFlowMassGridGeometry, couplingManager);
    using DarcyProblem = GetPropType<DarcyTypeTag, Properties::Problem>;
    auto darcyProblem = std::make_shared<DarcyProblem>(darcyGridGeometry, couplingManager);

    // the solution vector
    Traits::SolutionVector sol;
    sol[freeFlowMomentumIndex].resize(freeFlowMomentumGridGeometry->numDofs());
    sol[freeFlowMassIndex].resize(freeFlowMassGridGeometry->numDofs());
    sol[porousMediumIndex].resize(darcyGridGeometry->numDofs());
    freeFlowMassProblem->applyInitialSolution(sol[freeFlowMassIndex]);
    freeFlowMomentumProblem->applyInitialSolution(sol[freeFlowMomentumIndex]);
    darcyProblem->applyInitialSolution(sol[porousMediumIndex]);
    auto solOld = sol;

    // the grid variables
    using FreeFlowMomentumGridVariables = GetPropType<FreeFlowMomentumTypeTag, Properties::GridVariables>;
    auto freeFlowMomentumGridVariables = std::make_shared<FreeFlowMomentumGridVariables>(freeFlowMomentumProblem, freeFlowMomentumGridGeometry);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    using FreeFlowMassGridVariables = GetPropType<FreeFlowMassTypeTag, Properties::GridVariables>;
    auto freeFlowMassGridVariables = std::make_shared<FreeFlowMassGridVariables>(freeFlowMassProblem, freeFlowMassGridGeometry);
    using DarcyGridVariables = GetPropType<DarcyTypeTag, Properties::GridVariables>;
    auto darcyGridVariables = std::make_shared<DarcyGridVariables>(darcyProblem, darcyGridGeometry);

    couplingManager->init(freeFlowMomentumProblem, freeFlowMassProblem, darcyProblem,
                          std::make_tuple(freeFlowMomentumGridVariables, freeFlowMassGridVariables, darcyGridVariables),
                          sol, solOld);

    freeFlowMomentumGridVariables->init(sol[freeFlowMomentumIndex]);
    freeFlowMassGridVariables->init(sol[freeFlowMassIndex]);
    darcyGridVariables->init(sol[porousMediumIndex]);

    // We get some time loop parameters from the input file
    // and instantiate the time loop
    using Scalar = typename Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

        // the assembler for a stationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(freeFlowMomentumProblem, freeFlowMassProblem, darcyProblem),
                                                 std::make_tuple(freeFlowMomentumGridGeometry,
                                                                 freeFlowMassGridGeometry,
                                                                 darcyGridGeometry),
                                                 std::make_tuple(freeFlowMomentumGridVariables,
                                                                 freeFlowMassGridVariables,
                                                                 darcyGridVariables),
                                                 couplingManager, timeLoop, solOld);

    // the linear solver
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // initialize the vtk output module
    VtkOutputModule vtkWriterFF(*freeFlowMassGridVariables, sol[freeFlowMassIndex], freeFlowMassProblem->name());
    GetPropType<FreeFlowMassTypeTag, Properties::IOFields>::initOutputModule(vtkWriterFF); // Add model specific output fields
    vtkWriterFF.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<FreeFlowMassGridVariables>>());

    VtkOutputModule pmVtkWriter(*darcyGridVariables, sol[porousMediumIndex],  darcyProblem->name());
    GetPropType<DarcyTypeTag, Properties::IOFields>::initOutputModule(pmVtkWriter);
    pmVtkWriter.addVelocityOutput(std::make_shared<GetPropType<DarcyTypeTag, Properties::VelocityOutput>>(*darcyGridVariables));


    // write vtk output
    vtkWriterFF.write(0.0);
    pmVtkWriter.write(0.0);

        // timeloop
    timeLoop->start(); do
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(sol, *timeLoop);

        // make the new solution the old solution.
        solOld = sol;
        freeFlowMassGridVariables->advanceTimeStep();
        freeFlowMomentumGridVariables->advanceTimeStep();
        darcyGridVariables->advanceTimeStep();

        // advance to the time loop to the next step.
        timeLoop->advanceTimeStep();

        // write vtk output for each time step.
        vtkWriterFF.write(timeLoop->time());
        pmVtkWriter.write(timeLoop->time());

        // report statistics of this time step.
        timeLoop->reportTimeStep();

        // set a new dt as suggested by the newton solver for the next time step.
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

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
} // end main
