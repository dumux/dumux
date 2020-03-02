// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
// ## The main file
// This is the main file for the shallow water example. Here we can see the programme sequence and how the system is solved using newton's method.
// ### Includes
#include <config.h>

// Standard header file for C++, to get time and date information.
#include <ctime>

// Standard header file for C++, for in- and output.
#include <iostream>

// Dumux is based on DUNE, the Distributed and Unified Numerics Environment, which provides several grid managers and linear solvers. So we need some includes from that.
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>

// We need the following class to simplify the writing of dumux simulation data to VTK format.
#include <dumux/io/vtkoutputmodule.hh>
// In Dumux a property system is used to specify the model. For this, different properties are defined containing type definitions, values and methods. All properties are declared in the file properties.hh.
#include <dumux/common/properties.hh>
// The following file contains the parameter class, which manages the definition of input parameters by a default value, the inputfile or the command line.
#include <dumux/common/parameters.hh>
// The file dumuxmessage.hh contains the class defining the start and end message of the simulation.
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
// The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the different supported grid managers.
#include <dumux/io/grid/gridmanager.hh>
// We include the linear solver to be used to solve the linear system
#include <dumux/linear/amgbackend.hh>
// We include the nonlinear newtons method
#include <dumux/nonlinear/newtonsolver.hh>
// Further we include assembler, which assembles the linear systems for finite volume schemes (box-scheme, tpfa-approximation, mpfa-approximation)
#include <dumux/assembly/fvassembler.hh>
// We include the problem file which defines initial and boundary conditions to describe our example problem
#include "problem.hh"

// ### Beginning of the main function
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // We define the type tag for this problem
    using TypeTag = Properties::TTag::RoughChannel;

    // We initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // We print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // We parse command line arguments and input file
    Parameters::init(argc, argv);

    // ### Create the grid
    // A gridmanager tries to create the grid either from a grid file or the input file.
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // We compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // ### Setup and solving of the problem
    // #### Setup
    // We create and initialize the finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables.
    // We need the finite volume geometry to build up the subcontrolvolumes (scv) and subcontrolvolume faces (scvf) for each element of the grid partition.
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // In the problem, we define the boundary and initial conditions.
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // We initialize the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    // And then use the solutionvector to intialize the gridVariables.
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // We get some time loop parameters from the input file.
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // We intialize the vtk output module. Each model has a predefined model specific output with relevant parameters for that model.
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables,x, problem->name());
    // We add the analytical solution ("exactWaterDepth" and "exactVelocityX") to the predefined specific output.
    vtkWriter.addField(problem->getExactWaterDepth(), "exactWaterDepth");
    vtkWriter.addField(problem->getExactVelocityX(), "exactVelocityX");
    // We calculate the analytic solution.
    problem->analyticalSolution();
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    // We instantiate time loop.
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    //we set the assembler with the time loop because we have an instationary problem.
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop);

    // We set the linear solver.
    using LinearSolver = AMGBackend<TypeTag>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

    // Additionaly, we set the non-linear solver.
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // We set some check point at the end of the time loop. The check point is used to trigger the vtk output.
    timeLoop->setCheckPoint(tEnd);

    // We start the time loop.
    timeLoop->start(); do
    {
        // We start to calculate the new solution of that time step. First we define the old solution as the solution of the previous time step for storage evaluations.
        assembler->setPreviousSolution(xOld);

        // We solve the non-linear system with time step control.
        nonLinearSolver.solve(x,*timeLoop);

        // We make the new solution the old solution.
        xOld = x;
        gridVariables->advanceTimeStep();

        // We advance to the time loop to the next step.
        timeLoop->advanceTimeStep();

        // We write vtk output, if we reached the check point (end of time loop)
        if (timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());

        // We report statistics of this time step.
        timeLoop->reportTimeStep();

        // We set new dt as suggested by newton controller for the next time step.
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));


    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    // ### Final Output
    // We print dumux end message.
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main

catch (const Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (const Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (const Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
