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

 // ## The main file (`main.cc`)
 // [[content]]
 //
 // ### Included header files
 // [[details]] includes
 // [[exclude]]
 // Some generic includes.
#include <config.h>
#include <iostream>
// [[/exclude]]

// These is DUNE helper class related to parallel computation
#include <dune/common/parallel/mpihelper.hh>

// The following headers include functionality related to property definition or retrieval, as well as
// the retrieval of input parameters specified in the input file or via the command line.
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

// The following files contain the non-linear Newton solver, the available linear solver backends and the assembler for the linear
// systems arising from the staggered-grid discretization.
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/staggeredfvassembler.hh>

// The gridmanager constructs a grid from the information in the input or grid file.
// Many different Dune grid implementations are supported, of which a list can be found
// in `gridmanager.hh`.
#include <dumux/io/grid/gridmanager_yasp.hh>

// This class contains functionality for VTK output for models using the staggered finite volume scheme.
#include <dumux/io/staggeredvtkoutputmodule.hh>

// We include the problem header used for this simulation.
#include "properties.hh"
// [[/details]]
//
// The following function writes the velocities and coordinates at x = 0.5 and y = 0.5 into a log file.
// [[codeblock]]
template<class Problem, class SolutionVector, class GridGeometry>
void writeSteadyVelocityAndCoordinates(const Problem& problem, const SolutionVector &sol, const GridGeometry gridGeometry)
{
    std::ofstream logFilevx(problem->name() + "_vx.log"), logFilevy(problem->name() + "_vy.log");
    logFilevx << "y vx\n";
    logFilevy << "x vy\n";

    static constexpr double eps_ = 1.0e-7;
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        auto fvGeometry = localView(*gridGeometry);
        fvGeometry.bind(element);
        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (!scvf.boundary() && scvf.insideScvIdx() > scvf.outsideScvIdx())
            {
                const auto& globalPos = scvf.ipGlobal();
                const auto velocity = sol[gridGeometry->faceIdx()][scvf.dofIndex()][0];

                if (std::abs(globalPos[0]-0.5) < eps_)
                    logFilevx << globalPos[1] << " " << velocity << "\n";
                else if (std::abs(globalPos[1]-0.5) < eps_)
                    logFilevy << globalPos[0] << " " << velocity << "\n";
            }
        }
    }
}
// [[/codeblock]]

// ### The main function
// We will now discuss the main program flow implemented within the `main` function.
// At the beginning of each program using Dune, an instance of `Dune::MPIHelper` has to
// be created. Moreover, we parse the run-time arguments from the command line and the
// input file:
// [[codeblock]]
int main(int argc, char** argv)
{
    using namespace Dumux;

    // The Dune MPIHelper must be instantiated for each program using Dune, it is finalized automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);
    // [[/codeblock]]

    // We define a convenience alias for the type tag of the problem. The type
    // tag contains all the properties that are needed to define the model and the problem
    // setup. Throughout the main file, we will obtain types defined for this type tag
    // using the property system, i.e. with `GetPropType`.
    using TypeTag = Properties::TTag::LidDrivenCavityExample;

    // #### Step 1: Create the grid
    // The `GridManager` class creates the grid from information given in the input file.
    // This can either be a grid file, or in the case of structured grids, one can specify the coordinates
    // of the corners of the grid and the number of cells to be used to discretize each spatial direction.
    // [[codeblock]]
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // We compute on the leaf grid view.
    const auto& leafGridView = gridManager.grid().leafGridView();
    // [[/codeblock]]

    // #### Step 2: Setting up and solving the problem
    // First, a finite volume grid geometry is constructed from the grid that was created above.
    // This builds the sub-control volumes (scv) and sub-control volume faces (scvf) for each element
    // of the grid partition.
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // We now instantiate the problem, in which we define the boundary and initial conditions.
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // We set a solution vector which consist of two parts: one part (indexed by `cellCenterIdx`)
    // is for the pressure degrees of freedom (`dofs`) living in grid cell centers. Another part
    // (indexed by `faceIdx`) is for degrees of freedom defining the normal velocities on grid cell faces.
    // We initialize the solution vector by what was defined as the initial solution of the the problem.
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    x[GridGeometry::cellCenterIdx()].resize(gridGeometry->numCellCenterDofs());
    x[GridGeometry::faceIdx()].resize(gridGeometry->numFaceDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    // We use the initial solution vector to intialize the `gridVariables`.
    // The grid variables are used store variables (primary and secondary variables) on sub-control volumes and faces (volume and flux variables).
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // We get some time loop parameters from the input file
    // and instantiate the time loop
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // We then initialize the predefined model-specific output vtk output.
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    // To solve the non-linear problem at hand, we use the `NewtonSolver`,
    // which we have to tell how to assemble and solve the system in each
    // iteration. Here, we use the direct linear solver UMFPack.
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // ##### The time loop
    // In each time step, we solve the non-linear system of equations, write
    // the current solution into .vtk files and prepare for the next time step.
    // [[codeblock]]
    timeLoop->start(); do
    {
        // We solve the non-linear system with time step control.
        nonLinearSolver.solve(x, *timeLoop);

        // We make the new solution the old solution.
        xOld = x;
        gridVariables->advanceTimeStep();

        // We advance to the time loop to the next step.
        timeLoop->advanceTimeStep();

        // We write vtk output for each time step.
        vtkWriter.write(timeLoop->time());

        // We report statistics of this time step.
        timeLoop->reportTimeStep();

        // We set a new dt as suggested by the newton solver for the next time step.
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());
    // [[/codeblock]]

    // We write the velocities and coordinates at x = 0.5 and y = 0.5 into a file
    writeSteadyVelocityAndCoordinates(problem, x, gridGeometry);

    // The following piece of code prints a final status report of the time loop
    // before the program is terminated.
     // [[codeblock]]
    timeLoop->finalize(leafGridView.comm());

    // print used and unused parameters
    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;
} // end main
// [[/codeblock]]
// [[/content]]
