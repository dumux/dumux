// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TracerTests
 * \brief Test for the tracer model.
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/pdesolver.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/experimental/assembly/multistagefvassembler.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    //! define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // parse the command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // setup instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    //! we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    //! create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    //! the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    //! the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    //! the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    //! initialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVolumeVariable([](const auto& v){ return v.saturation(); }, "S_aq");
    vtkWriter.write(0.0);

    //! get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = getParam<Scalar>("TimeLoop.Dt");

    //! instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    auto timeSteppingMethod = std::make_shared<Experimental::MultiStage::ExplicitEuler<Scalar>>();

    //! the assembler with time loop for instationary problem
    using Assembler = Experimental::MultiStageFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeSteppingMethod, xOld);

    //! the linear solver
    using LinearSolver = ExplicitDiagonalSolver;
    auto linearSolver = std::make_shared<LinearSolver>();

    //! the system solver
    using Solver = LinearPDESolver<Assembler, LinearSolver>;
    auto solver = std::make_shared<Solver>(assembler, linearSolver);

    using TimeStepper = Experimental::MultiStageTimeStepper<Solver>;
    TimeStepper timeStepper(solver, timeSteppingMethod);

    //! set some check points for the time loop
    timeLoop->setPeriodicCheckPoint(tEnd/10.0);

    //! start the time loop
    timeLoop->start();
    while (!timeLoop->finished())
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

        // write vtk output on check points
        if (timeLoop->isCheckPoint() || timeLoop->finished())
            vtkWriter.write(timeLoop->time());

        // set new dt TODO implement CFL-criterion
        timeLoop->setTimeStepSize(dt);
    }

    timeLoop->finalize(leafGridView.comm());

    if (leafGridView.comm().rank() == 0)
        Parameters::print();

    return 0;

}
