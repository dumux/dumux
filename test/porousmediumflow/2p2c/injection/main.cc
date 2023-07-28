// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTwoCTests
 * \brief Test for the problem where air is injected under a low permeable layer in a depth of 2700m.
 */
#include <config.h>

#include <iostream>

#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/loadsolution.hh>

#include <dumux/experimental/assembly/multistagefvassembler.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    if (restartTime > 0)
    {
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        const auto fileName = getParam<std::string>("Restart.File");
        const auto pvName = createPVNameFunction<IOFields, PrimaryVariables, ModelTraits, FluidSystem>();
        loadSolution(x, fileName, pvName, *gridGeometry);
    }
    else
        problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // initialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(restartTime);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    auto timeSteppingMethod = std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>();

    // the assembler with time loop for instationary problem
    using Assembler = Experimental::MultiStageFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeSteppingMethod, xOld);

    // the linear solver
    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver);

    using TimeStepper = Experimental::MultiStageTimeStepper<NewtonSolver>;
    TimeStepper timeStepper(nonLinearSolver, timeSteppingMethod);

    // time loop
    timeLoop->start(); do
    {
        // time integration
        timeStepper.step(x, timeLoop->time(), timeLoop->timeStepSize());

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize()));

        // write vtk output
        vtkWriter.write(timeLoop->time());

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    if (leafGridView.comm().rank() == 0)
        Parameters::print();

    return 0;
} // end main
