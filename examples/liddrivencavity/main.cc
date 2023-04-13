// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

 // ## The main file (`main.cc`)
 // [[content]]
 //
 // ### Included header files
 // [[details]] includes
 // [[exclude]]
 // Some generic includes.
#include <config.h>
#include <iostream>
#include <dune/common/timer.hh>
#include <dumux/common/dumuxmessage.hh>
// [[/exclude]]

// These is DUNE helper class related to parallel computation
#include <dune/common/parallel/mpihelper.hh>

// The following headers include functionality related to property definition or retrieval, as well as
// the retrieval of input parameters specified in the input file or via the command line.
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/initialize.hh>

// The following files contain the multi-domain Newton solver, the available linear solver backends and the assembler for the linear
// systems arising from the staggered-grid discretization.
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

// The gridmanager constructs a grid from the information in the input or grid file.
// Many different Dune grid implementations are supported, of which a list can be found
// in `gridmanager.hh`.
#include <dumux/io/grid/gridmanager_yasp.hh>

// This class contains functionality for VTK output for models using the staggered finite volume scheme.
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>

// We include the problem header used for this simulation.
#include "properties.hh"
// [[/details]]

// The following function writes the velocities and coordinates at x = 0.5 and y = 0.5 into a log file.
// [[codeblock]]
template<class Problem, class SolutionVector>
void writeSteadyVelocityAndCoordinates(const Problem& problem, const SolutionVector& sol)
{
    const auto& gridGeometry = problem.gridGeometry();
    std::ofstream logFilevx(problem.name() + "_vx.log"), logFilevy(problem.name() + "_vy.log");
    logFilevx << "y vx\n";
    logFilevy << "x vy\n";

    static constexpr double eps_ = 1.0e-7;
    std::vector<int> dofHandled(gridGeometry.numDofs(), false);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        auto fvGeometry = localView(gridGeometry);
        fvGeometry.bind(element);
        for (const auto& scv : scvs(fvGeometry))
        {
            if (dofHandled[scv.dofIndex()])
                continue;

            if (!scv.boundary())
            {
                const auto& globalPos = scv.dofPosition();
                const auto velocity = sol[scv.dofIndex()][0];
                dofHandled[scv.dofIndex()] = true;

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

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // parse command line arguments and input file
    Parameters::init(argc, argv);
    // [[/codeblock]]

    // We define a convenience alias for the type tags of the problems. The type
    // tag contains all the properties that are needed to define the model and the problem
    // setup. Throughout the main file, we will obtain types defined each type tag
    // using the property system, i.e. with `GetPropType`.
    using MomentumTypeTag = Properties::TTag::LidDrivenCavityExampleMomentum;
    using MassTypeTag = Properties::TTag::LidDrivenCavityExampleMass;

    // #### Step 1: Create the grid
    // The `GridManager` class creates the grid from information given in the input file.
    // This can either be a grid file, or in the case of structured grids, one can specify the coordinates
    // of the corners of the grid and the number of cells to be used to discretize each spatial direction.
    // [[codeblock]]
    GridManager<GetPropType<MomentumTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // We compute on the leaf grid view.
    const auto& leafGridView = gridManager.grid().leafGridView();
    // [[/codeblock]]

    // #### Step 2: Setting up and solving the problem
    // First, a finite volume grid geometry is constructed from the grid that was created above.
    // This builds the sub-control volumes (scv) and sub-control volume faces (scvf) for each element
    // of the grid partition.
    // This is done for both the momentum and mass grid geometries
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // We introduce the multidomain coupling manager, which will couple the mass and the momentum problems
    // We can obtain the type from either the `MomentumTypeTag` or the `MassTypeTag` because they are mutually coupled with the same manager
    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // We now instantiate the problems, in which we define the boundary and initial conditions.
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // We set a solution vector which consist of two parts: one part (indexed by `massIdx`)
    // is for the pressure degrees of freedom (`dofs`) living in grid cell centers. Another part
    // (indexed by `momentumIdx`) is for degrees of freedom defining the normal velocities on grid cell faces.
    // We initialize the solution vector by what was defined as the initial solution of the the problem.
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);
    auto xOld = x;

    // We use the initial solution vector to create the `gridVariables`.
    // The grid variables are used store variables (primary and secondary variables) on sub-control volumes and faces (volume and flux variables).
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    // using the problems and the grid variables, the coupling manager and the grid variables are initialized with the initial solution.
    // The grid variables have to be initialized _after_ the coupling manager.
    // This is because they require the correct coupling context between the mass and momentum model to be initialized.
    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);
    momentumGridVariables->init(x[momentumIdx]);
    massGridVariables->init(x[massIdx]);

    // We get some time loop parameters from the input file
    // and instantiate the time loop
    using Scalar = typename Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // We then initialize the predefined model-specific output vtk output.
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    vtkWriter.write(0.0);

    // To solve the non-linear problem at hand, we use the `NewtonSolver`,
    // which we have to tell how to assemble and solve the system in each
    // iteration. Here, we use the direct linear solver UMFPack.
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager, timeLoop, xOld);
    // the linear solver
    using LinearSolver = Dumux::UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

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
        massGridVariables->advanceTimeStep();
        momentumGridVariables->advanceTimeStep();

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
    writeSteadyVelocityAndCoordinates(*momentumProblem, x[momentumIdx]);

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
