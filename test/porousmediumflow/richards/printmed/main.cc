// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsTests
 * \brief Test for the Richards model.
 */

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/linear/istlsolverfactorybackend.hh>

#include <dumux/porousmediumflow/richards/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/loadsolution.hh>

#include "properties.hh"
#include "dropsolver.hh"

#ifndef DIFFMETHOD
#define DIFFMETHOD DiffMethod::numeric
#endif

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::RichardsLensBox;

    // initialize Dumux (parallel helpers)
    // always call this before any other code
    Dumux::initialize(argc, argv);

    // get an instance of the MPI helper
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());

    using DropSolver = DropletSolverTwoP<TypeTag, false>;
    auto dropSolver = std::make_shared<DropSolver>(*problem, x);
    problem->setDropSolver(dropSolver);

    problem->applyInitialSolution(x);
    auto xOld = x;

    // // maybe update the interface parameters
    // if (ENABLEINTERFACESOLVER)
    //     problem->spatialParams().updateMaterialInterfaces(x);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);
    dropSolver->setGridVariables(gridVariables);

    // initialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;

    // use non-conforming output for the test with interface solver
    const auto ncOutput = getParam<bool>("Problem.UseNonConformingOutput", false);
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name(), "",
                                                             ncOutput ? Dune::VTK::nonconforming : Dune::VTK::conforming);
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(restartTime);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    dropSolver->setTimeLoop(timeLoop);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DIFFMETHOD>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = IstlSolverFactoryBackend<LinearSolverTraits<GridGeometry>,
                                                  LinearAlgebraTraitsFromAssembler<Assembler>>;

    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::RichardsNewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        dropSolver->update();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    // write Newton report
    nonLinearSolver.report();

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