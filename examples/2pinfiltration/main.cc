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
// ## The file `main.cc`
// [[content]]
//
// This is the main file for the 2pinfiltration example. Here we can see the programme sequence and how the system is solved using Newton's method
// ### Included header files
// [[details]] includes
// [[codeblock]]
#include <config.h>

// Standard header file for C++, for in- and output.
#include <iostream>
// [[/codeblock]]

// A Dune helper class for enabling parallel simulations with MPI
#include <dune/common/parallel/mpihelper.hh>

// In Dumux, a property system is used to specify the model. For this, different properties are defined containing type definitions, values and methods. All properties are declared in the file properties.hh. Additionally, we include the parameter class, which manages the definition of input parameters by a default value, the inputfile or the command line.
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

//We include the linear solver to be used to solve the linear system and the nonlinear  Newton's method
#include <dumux/linear/amgbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

// Further, we include assembler, which assembles the linear systems for finite volume schemes (box-scheme, tpfa-approximation, mpfa-approximation) and a file that defines the different differentiation methods used to compute the derivatives of the residual
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

// We need the following class to simplify the writing of dumux simulation data to VTK format.
#include <dumux/io/vtkoutputmodule.hh>
// The gridmanager constructs a grid from the information in the input or grid file.
#include <dumux/io/grid/gridmanager_alu.hh>

//We include several files which are needed for the adaptive grid
#include <dumux/adaptive/adapt.hh>
#include <dumux/adaptive/markelements.hh>
#include <dumux/adaptive/initializationindicator.hh>
#include <dumux/porousmediumflow/2p/griddatatransfer.hh>
#include <dumux/porousmediumflow/2p/gridadaptindicator.hh>

// Finally, we include the properties which configure the simulation
#include "properties.hh"
// [[/details]]
//

// ### The main function
// At the beginning of the simulation, we create a type tag with a suitable alias for our problem. // This will contain all the properties that are needed to define the problem set-up and the model we use. Additionally, we have to create an instance of `Dune::MPIHelper` and parse the run-time arguments:
// [[codeblock]]
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // we define the type tag for this problem
    using TypeTag = Properties::TTag::PointSourceExample;

    //We initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    //We parse command line arguments and input file
    Parameters::init(argc, argv);
    // [[/codeblock]]

    // #### Create the grid
    // A gridmanager tries to create the grid either from a grid file or the input file.
    // [[codeblock]]
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // The instationary non-linear problem is run on this grid.
    //
    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();
    // [[/codeblock]]

    // #### Set-up of the problem
    //
    // We create and initialize the finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables.
    //
    // We need the finite volume geometry to build up the subcontrolvolumes (scv) and subcontrolvolume faces (scvf) for each element of the grid partition.
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // In the problem, we define the boundary and initial conditions and compute the point sources. The `computePointSourceMap` method is inherited from the fvproblem and therefore specified in the `dumux/common/fvproblem.hh`. It calls the `addPointSources` method specified in the `problem.hh` file.
     // [[codeblock]]
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);
    // We call the `computePointSourceMap` method to compute the point sources.
    problem->computePointSourceMap();
     // [[/codeblock]]

    // We initialize the solution vector and then use the solution vector to intialize the `gridVariables`.
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // ##### Grid adaption
    // We instantiate the indicator for grid adaption & the data transfer, we read some parameters for indicator from the input file.
    // [[codeblock]]
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const Scalar refineTol = getParam<Scalar>("Adaptive.RefineTolerance");
    const Scalar coarsenTol = getParam<Scalar>("Adaptive.CoarsenTolerance");
    // We use an indicator for a two-phase flow problem that is saturation-dependent.
    // It is defined in the file `dumux/porousmediumflow/2p/gridadaptindicator.hh.`
    // and allows to set the minimum and maximum allowed refinement levels via the input parameters.
    TwoPGridAdaptIndicator<TypeTag> indicator(gridGeometry);
    // The data transfer performs the transfer of data on a grid from before to after adaptation
    // and is defined in the file `dumux/porousmediumflow/2p/griddatatransfer.hh`.
    // Its main functions are to store and reconstruct the primary variables.
    TwoPGridDataTransfer<TypeTag> dataTransfer(problem, gridGeometry, gridVariables, x);
    // [[/codeblock]]

    // We do an initial refinement around sources/BCs. We use the `GridAdaptInitializationIndicator` defined in `dumux/adaptive/initializationindicator.hh` for that.
    GridAdaptInitializationIndicator<TypeTag> initIndicator(problem, gridGeometry, gridVariables);

    //We refine up to the maximum level. For every level, the indicator used for the refinement/coarsening is calculated. If any grid cells have to be adapted, the gridvariables and the pointsourcemap are updated.
     // [[codeblock]]
    const auto maxLevel = getParam<std::size_t>("Adaptive.MaxLevel", 0);
    for (std::size_t i = 0; i < maxLevel; ++i)
    {
        //we calculate the initial indicator for adaption for each grid cell using the initial solution x
        initIndicator.calculate(x);

        //and then we mark the elements that were adapted.
        bool wasAdapted = false;
        if (markElements(gridManager.grid(), initIndicator))
            wasAdapted = adapt(gridManager.grid(), dataTransfer);

        // In case of a grid adaptation, the gridvariables and the pointsourcemap are updated.
        if (wasAdapted)
        {
            // We overwrite the old solution with the new (resized & interpolated) one
            xOld = x;
            //We initialize the secondary variables to the new (and "new old") solution
            gridVariables->updateAfterGridAdaption(x);
            // we update the point source map after adaption
            problem->computePointSourceMap();
        }
    }
    // [[/codeblock]]

    // Depending on the initial conditions, another grid adaptation might be necessary.
    //The gridadaptindicator uses the input parameters `Adaptive.RefineTolerance` and `Adaptive.CoarsenTolerance` for this step.
    //Again, if elements were adapted, we mark them and update gridvariables and the pointsourcemap accordingly.
    // [[codeblock]]
    indicator.calculate(x, refineTol, coarsenTol);

    //we mark the elements that were adapted
    bool wasAdapted = false;
    if (markElements(gridManager.grid(), indicator))
        wasAdapted = adapt(gridManager.grid(), dataTransfer);

    // In case of an additional grid adaptation, the gridvariables and the pointsourcemap are updated again.
    if (wasAdapted)
    {
        // Overwrite the old solution with the new (resized & interpolated) one
        xOld = x;
        // Initialize the secondary variables to the new (and "new old") solution
        gridVariables->updateAfterGridAdaption(x);
        // Update the point source map
        problem->computePointSourceMap();
    }
    // [[/codeblock]]

    // #### Solving the problem

    // We get some time loop parameters from the input file params.input
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // and initialize the vtkoutput. Each model has a predefined model specific output with relevant parameters for that model.
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    // We instantiate the time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // and set the assembler with the time loop because we have an instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // We set the linear solver and the non-linear solver
    using LinearSolver = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // ##### The time loop
    // We start the time loop. In each time step before we start calculating a new solution we check if we have to refine the mesh again based on the solution of the previous time step. If the grid is adapted we update the gridvariables and the pointsourcemap. Afterwards, the solution of that time step is calculated.
    // [[codeblock]]
    timeLoop->start(); do
    {
        // We only want to refine/coarsen after first time step is finished, not before.
        //The initial refinement was already done before the start of the time loop.
        //This means we only refine when the time is greater than 0.
        if (timeLoop->time() > 0)
        {
            // again we compute the refinement indicator with the `TwoPGridAdaptIndicator`
            indicator.calculate(x, refineTol, coarsenTol);

            //we mark elements and adapt grid if necessary
            wasAdapted = false;
            if (markElements(gridManager.grid(), indicator))
                wasAdapted = adapt(gridManager.grid(), dataTransfer);

            //In case of a grid adaptation, the gridvariables and the pointsourcemap are updated again.
            if (wasAdapted)
            {
                // We overwrite the old solution with the new (resized & interpolated) one
                xOld = x;
                 // We tell the assembler to resize the matrix and set pattern
                assembler->setJacobianPattern();
                // We tell the assembler to resize the residual
                assembler->setResidualSize();
                 // We initialize the secondary variables to the new (and "new old") solution
                gridVariables->updateAfterGridAdaption(x);
                // We update the point source map
                problem->computePointSourceMap();
            }
        // we leave the refinement step
        }

        // We solve the non-linear system with time step control.
        nonLinearSolver.solve(x, *timeLoop);

        //We make the new solution the old solution.
        xOld = x;
        gridVariables->advanceTimeStep();

        //We advance to the time loop to the next step.
        timeLoop->advanceTimeStep();

        //We write vtk output for each time step
        vtkWriter.write(timeLoop->time());

        //We report statistics of this time step
        timeLoop->reportTimeStep();

        //We set a new dt as suggested by the newton solver for the next time step
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());
    // [[/codeblock]]

    // The following piece of code prints a final status report of the time loop
    //  before the program is terminated and we print he dumux end message
     // [[codeblock]]
    timeLoop->finalize(leafGridView.comm());

    // print dumux end message
    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;
} // end main
//[[/codeblock]]
// ### Exception handling
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
// [[/codeblock]]
// [[/details]]
// [[/content]]
