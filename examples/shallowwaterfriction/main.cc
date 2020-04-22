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

// ## The main file `main.cc`
// [[content]]
//
// ### Included header files
// [[details]] includes
// [[exclude]]
// Some generic includes.
#include <config.h>
#include <iostream>
// [[/exclude]]
//
// DUNE helper class for MPI
#include <dune/common/parallel/mpihelper.hh>

// The following headers include functionality related to property definition or retrieval, as well as
// the retrieval of input parameters specified in the input file or via the command line.
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

// The following files contains the available linear solver backends, the non linear Newton Solver
// and the assembler for the linear systems arising from finite volume discretizations
// (box-scheme, tpfa-approximation, mpfa-approximation).
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

// The following class provides a convenient way of writing of dumux simulation results to VTK format.
#include <dumux/io/vtkoutputmodule.hh>
// The gridmanager constructs a grid from the information in the input or grid file.
// Many different Dune grid implementations are supported, of which a list can be found
// in `gridmanager.hh`.
#include <dumux/io/grid/gridmanager_yasp.hh>

// We include the header file specifing the properties of this example
#include "properties.hh"
// [[/details]]
//

// ### The main function
// We will now discuss the main program flow implemented within the `main` function.
// At the beginning of each program using Dune, an instance of `Dune::MPIHelper` has to
// be created. Moreover, we parse the run-time arguments from the command line and the
// input file:
// [[codeblock]]
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // The Dune MPIHelper must be instantiated for each program using Dune
    Dune::MPIHelper::instance(argc, argv);

    // We parse command line arguments and input file
    Parameters::init(argc, argv);
    // [[/codeblock]]

    // We define a convenience alias for the type tag of the problem. The type
    // tag contains all the properties that are needed to define the model and the problem
    // setup. Throughout the main file, we will obtain types defined for these type tag
    // using the property system, i.e. with `GetPropType`.
    using TypeTag = Properties::TTag::RoughChannel;

    // #### Step 1: Create the grid
    // The `GridManager` class creates the grid from information given in the input file.
    // This can either be a grid file or, in the case of structured grids, one can specify the coordinates
    // of the corners of the grid and the number of cells to be used to discretize each spatial direction.
    //[[codeblock]]
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // We compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    //[[/codeblock]]

    // #### Step 2: Solving the shallow water problem
    // First, a finite volume grid geometry is constructed from the grid that was created above.
    // This builds the sub-control volumes (scv) and sub-control volume faces (scvf) for each element
    // of the grid partition.
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // We now instantiate the problem, in which we define the boundary and initial conditions.
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // We initialize the solution vector. The shallow water problem is transient,
    // therefore the initial solution defined in the problem is applied to the
    // solution vector. On the basis of this solution, we initialize then the grid variables.
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // Let us now instantiate the time loop. Therefore, we read in some time loop parameters from the input file.
    // The parameter `tEnd` defines the duration of the simulation, `dt` the initial time step size and `maxDt` the maximal time step size.
    // Moreover, we define the end of the simulation `tEnd` as check point in the time loop at which we will write the solution to vtk files.
    // [[codeblock]]
    using Scalar = GetPropType<TypeTag, Properties::Scalar>; // type for scalar values
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // We instantiate time loop.
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    timeLoop->setCheckPoint(tEnd);
    // [[/codeblock]]

    // We initialize the assembler with a time loop for the transient problem.
    // Within the time loop, we will use this assembler in each time step to assemble the linear system.
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // We initialize the linear solver.
    using LinearSolver = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

    // We initialize the non-linear solver.
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // The following lines of code initialize the vtk output module, add the velocity output facility
    // and write out the initial solution. At each checkpoint, we will use the output module to write
    // the solution of a time step into a corresponding vtk file.
    // [[codeblock]]
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables,x, problem->name());

    // add model-specific output fields to the writer
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);

    // We add the analytical solution ("exactWaterDepth" and "exactVelocityX") to the predefined specific output.
    vtkWriter.addField(problem->getExactWaterDepth(), "exactWaterDepth");
    vtkWriter.addField(problem->getExactVelocityX(), "exactVelocityX");

    // We calculate the analytic solution.
    problem->analyticalSolution();

    // write initial solution (including the above calculated analytical solution.
    vtkWriter.write(0.0);
    // [[/codeblock]]

    // ##### The time loop
    // We start the time loop and solve a new time step as long as `tEnd` is not reached. In every time step,
    // the problem is assembled and solved, the solution is updated, and when a checkpoint is reached the solution
    // is written to a new vtk file. In addition, statistics related to CPU time, the current simulation time
    // and the time step sizes used is printed to the terminal.
    // [[codeblock]]
    timeLoop->start(); do
    {
        // We solve the non-linear system with time step control, using Newthon's method.
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
    // [[/codeblock]]

    // The following piece of code prints a final status report of the time loop
    //  before the program is terminated.
    timeLoop->finalize(leafGridView.comm());

    return 0;
}
// #### Exception handling
// In this part of the main file we catch and print possible exceptions that could
// occur during the simulation.
// [[details]] exception handler
// [[codeblock]]
// errors related to run-time parameters
catch (const Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
// errors related to the parsing of Dune grid files
catch (const Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
// generic error handling with Dune::Exception
catch (const Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
// other exceptions
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
// [[/codeblock]]
// [[/details]]
// [[/content]]
