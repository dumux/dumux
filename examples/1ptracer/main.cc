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

// We look now at the main file for the tracer problem. We set up two problems in this file and solve them sequentially, first the 1p problem and afterwards the tracer problem. The result of the 1p problem is the pressure distribution in the problem domain. We use it to calculate the volume fluxes, which act as an input for the tracer problem. Based on this volume fluxes, we calculate the transport of a tracer in the following tracer problem.
// ### Includes
#include <config.h>

// We include both problems in the main file, the `problem_1p.hh` and the `problem_tracer.hh`.
#include "problem_1p.hh"
#include "problem_tracer.hh"

// Further, we include a standard header file for C++, to get time and date information
#include <ctime>
// and another one for in- and output.
#include <iostream>

// Dumux is based on DUNE, the Distributed and Unified Numerics Environment, which provides several grid managers and linear solvers. So we need some includes from that.
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>

// In Dumux, a property system is used to specify the model. For this, different properties are defined containing type definitions, values and methods. All properties are declared in the file `properties.hh`.
#include <dumux/common/properties.hh>
// The following file contains the parameter class, which manages the definition of input parameters by a default value, the inputfile or the command line.
#include <dumux/common/parameters.hh>
// The file `dumuxmessage.hh` contains the class defining the start and end message of the simulation.
#include <dumux/common/dumuxmessage.hh>

// The following file contains the class, which defines the sequential linear solver backends.
#include <dumux/linear/seqsolverbackend.hh>
// Further we include the assembler, which assembles the linear systems for finite volume schemes (box-scheme, tpfa-approximation, mpfa-approximation).
#include <dumux/assembly/fvassembler.hh>
// The containing class in the following file defines the different differentiation methods used to compute the derivatives of the residual.
#include <dumux/assembly/diffmethod.hh>

// We need the following class to simplify the writing of dumux simulation data to VTK format.
#include <dumux/io/vtkoutputmodule.hh>
// The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the different supported grid managers.
#include <dumux/io/grid/gridmanager.hh>

// ### Beginning of the main function
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // We define the type tags for the two problems. They are created in the individual problem files.
    using OnePTypeTag = Properties::TTag::IncompressibleTest;
    using TracerTypeTag = Properties::TTag::TracerTestCC;

    // We initialize MPI. Finalization is done automatically on exit.
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // We print the dumux start message.
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // We parse the command line arguments.
    Parameters::init(argc, argv);

    // ### Create the grid
    // A gridmanager tries to create the grid either from a grid file or the input file. We solve both problems on the same grid. Hence, the grid is only created once using the type tag of the 1p problem.
    GridManager<GetPropType<OnePTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // We compute on the leaf grid view.
    const auto& leafGridView = gridManager.grid().leafGridView();

    // ### Set-up and solving of the 1p problem
    // In the following section, we set up and solve the 1p problem. As the result of this problem, we obtain the pressure distribution in the problem domain.
    // #### Set-up
    // We create and initialize the finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables.
    // We need the finite volume geometry to build up the subcontrolvolumes (scv) and subcontrolvolume faces (scvf) for each element of the grid partition.
    using GridGeometry = GetPropType<OnePTypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // In the problem, we define the boundary and initial conditions.
    using OnePProblem = GetPropType<OnePTypeTag, Properties::Problem>;
    auto problemOneP = std::make_shared<OnePProblem>(gridGeometry);

    // The jacobian matrix (`A`), the solution vector (`p`) and the residual (`r`) are parts of the linear system.
    using JacobianMatrix = GetPropType<OnePTypeTag, Properties::JacobianMatrix>;
    using SolutionVector = GetPropType<OnePTypeTag, Properties::SolutionVector>;
    SolutionVector p(leafGridView.size(0));

    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<SolutionVector>();

    // The grid variables store variables on scv and scvf (volume and flux variables).
    using OnePGridVariables = GetPropType<OnePTypeTag, Properties::GridVariables>;
    auto onePGridVariables = std::make_shared<OnePGridVariables>(problemOneP, gridGeometry);
    onePGridVariables->init(p);

    // #### Assembling the linear system
    // We created and inizialize the assembler.
    using OnePAssembler = FVAssembler<OnePTypeTag, DiffMethod::analytic>;
    auto assemblerOneP = std::make_shared<OnePAssembler>(problemOneP, gridGeometry, onePGridVariables);
    assemblerOneP->setLinearSystem(A, r);

    // We assemble the local jacobian and the residual and stop the time needed, which is displayed in the terminal output, using the `assemblyTimer`. Further, we start the timer to evaluate the total time of the assembly, solving and updating.
    Dune::Timer timer;
    Dune::Timer assemblyTimer; std::cout << "Assembling linear system ..." << std::flush;
    assemblerOneP->assembleJacobianAndResidual(p);
    assemblyTimer.stop(); std::cout << " took " << assemblyTimer.elapsed() << " seconds." << std::endl;

    // We want to solve `Ax = -r`.
    (*r) *= -1.0;

    // #### Solution
    // We set the linear solver "UMFPack" as the linear solver. Afterwards we solve the linear system. The time needed to solve the system is recorded by the `solverTimer` and displayed in the terminal output.
    using LinearSolver = UMFPackBackend;
    Dune::Timer solverTimer; std::cout << "Solving linear system ..." << std::flush;
    auto linearSolver = std::make_shared<LinearSolver>();
    linearSolver->solve(*A, p, *r);
    solverTimer.stop(); std::cout << " took " << solverTimer.elapsed() << " seconds." << std::endl;

    // #### Update and output
    // We update the grid variables with the new solution.
    Dune::Timer updateTimer; std::cout << "Updating variables ..." << std::flush;
    onePGridVariables->update(p);
    updateTimer.elapsed(); std::cout << " took " << updateTimer.elapsed() << std::endl;


    // We initialize the vtkoutput. Each model has a predefined model specific output with relevant parameters for that model. We add the pressure data from the solution vector (`p`) and the permeability field as output data.
    using GridView = GetPropType<OnePTypeTag, Properties::GridView>;
    Dune::VTKWriter<GridView> onepWriter(leafGridView);
    onepWriter.addCellData(p, "p");
    const auto& k = problemOneP->spatialParams().getKField();
    onepWriter.addCellData(k, "permeability");
    onepWriter.write("1p");

    // We stop the timer and display the total time of the simulation as well as the cumulative CPU time.
    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";


    // ### Computation of the volume fluxes
    // We use the results of the 1p problem to calculate the volume fluxes in the model domain.

    using Scalar =  GetPropType<OnePTypeTag, Properties::Scalar>;
    std::vector<Scalar> volumeFlux(gridGeometry->numScvf(), 0.0);

    using FluxVariables =  GetPropType<OnePTypeTag, Properties::FluxVariables>;
    auto upwindTerm = [](const auto& volVars) { return volVars.mobility(0); };

    // We iterate over all elements.
    for (const auto& element : elements(leafGridView))
    {
        auto fvGeometry = localView(*gridGeometry);
        fvGeometry.bind(element);

        auto elemVolVars = localView(onePGridVariables->curGridVolVars());
        elemVolVars.bind(element, fvGeometry, p);

        auto elemFluxVars = localView(onePGridVariables->gridFluxVarsCache());
        elemFluxVars.bind(element, fvGeometry, elemVolVars);

        // We calculate the volume flux for every subcontrolvolume face, which is not on a Neumann boundary (is not on the boundary or is on a Dirichlet boundary).

        for (const auto& scvf : scvfs(fvGeometry))
        {
            const auto idx = scvf.index();

            if (!scvf.boundary())
            {
                FluxVariables fluxVars;
                fluxVars.init(*problemOneP, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                volumeFlux[idx] = fluxVars.advectiveFlux(0, upwindTerm);
            }
            else
            {
                const auto bcTypes = problemOneP->boundaryTypes(element, scvf);
                if (bcTypes.hasOnlyDirichlet())
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*problemOneP, element, fvGeometry, elemVolVars, scvf, elemFluxVars);
                    volumeFlux[idx] = fluxVars.advectiveFlux(0, upwindTerm);
                }
            }
        }
    }


    // ### Set-up and solving of the tracer problem
    // #### Set-up
    // Similar to the 1p problem, we first create and initialize the problem.
    using TracerProblem = GetPropType<TracerTypeTag, Properties::Problem>;
    auto tracerProblem = std::make_shared<TracerProblem>(gridGeometry);

    // We use the volume fluxes calculated in the previous section as input for the tracer model.
    tracerProblem->spatialParams().setVolumeFlux(volumeFlux);

    // We create and initialize the solution vector. As the tracer problem is transient, the initial solution defined in the problem is applied to the solution vector.
    SolutionVector x(leafGridView.size(0));
    tracerProblem->applyInitialSolution(x);
    auto xOld = x;

    // We create and initialize the grid variables.
    using GridVariables = GetPropType<TracerTypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(tracerProblem, gridGeometry);
    gridVariables->init(x);

    // We read in some time loop parameters from the input file. The parameter `tEnd` defines the duration of the simulation, dt the initial time step size and `maxDt` the maximal time step size.
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");

    // We instantiate the time loop.
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // We create and inizialize the assembler with time loop for the instationary problem.
    using TracerAssembler = FVAssembler<TracerTypeTag, DiffMethod::analytic, /*implicit=*/false>;
    auto assembler = std::make_shared<TracerAssembler>(tracerProblem, gridGeometry, gridVariables, timeLoop);
    assembler->setLinearSystem(A, r);

    // We initialize the vtk output module and add a velocity output.
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, tracerProblem->name());
    using IOFields = GetPropType<TracerTypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    using VelocityOutput = GetPropType<TracerTypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    vtkWriter.write(0.0);


    // For the time loop we set a check point.
    timeLoop->setPeriodicCheckPoint(tEnd/10.0);

    // #### The time loop
    // We start the time loop and calculate a new time step as long as `tEnd` is not reached. In every single time step, the problem is assembled and solved.
    timeLoop->start(); do
    {
        // First we define the old solution as the solution of the previous time step for storage evaluations.
        assembler->setPreviousSolution(xOld);

        // Then the linear system is assembled.
        Dune::Timer assembleTimer;
        assembler->assembleJacobianAndResidual(x);
        assembleTimer.stop();

        // We solve the linear system `A(xOld-xNew) = r`.
        Dune::Timer solveTimer;
        SolutionVector xDelta(x);
        linearSolver->solve(*A, xDelta, *r);
        solveTimer.stop();

        // We calculate the actual solution and update it in the grid variables.
        updateTimer.reset();
        x -= xDelta;
        gridVariables->update(x);
        updateTimer.stop();

        // We display the statistics of the actual time step.
        const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
        std::cout << "Assemble/solve/update time: "
                  <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                  <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                  <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                  <<  std::endl;

        // The new solution is defined as the old solution.
        xOld = x;
        gridVariables->advanceTimeStep();

        // We advance the time loop to the next time step.
        timeLoop->advanceTimeStep();

        // We write the Vtk output on check points.
        if (timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());

        // We report the statistics of this time step.
        timeLoop->reportTimeStep();

        // We set the time step size dt of the next time step.
        timeLoop->setTimeStepSize(dt);

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());


    // ### Final Output

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;

}
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
