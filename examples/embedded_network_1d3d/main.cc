// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// # The main file (`main.cc`)
//
// This file contains the `main` function and implements the main program
// flow. We initialize the simulation framework, initialize parameters,
// create the grids (using parameters from the configuration file `params.input`),
// create vector to store the primary and secondary variables, construct an
// assembler that can assembly the the residual (discrete equations) and
// the system matrix (Jacobian of the residual), and create a linear solver
// that solves the linear system arising in each time step. The time loop
// implements a simple time stepping scheme with constant time step size.
//
// [[content]]
//
// ### Included header files
// [[details]] includes
// [[codeblock]]
#include <config.h>

#include <vector>
#include <numeric>

#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/initialize.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/container.hh>
#include <dumux/io/format.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include "bloodflow.hh"
#include "solver.hh"
#include "properties.hh"
// [[/codeblock]]
// [[/details]]
//
// ## The main function
// We will now discuss the main program flow implemented within the `main` function.
// At the beginning of each program using Dune, an instance of `Dune::MPIHelper` has to
// be created. Moreover, we parse the run-time arguments from the command line and the
// input file:
int main(int argc, char** argv) {
    // enable all symbols of the Dumux namespace without having to prefix with ``Dumux::``
    using namespace Dumux;
    //
    // ### initialization and parameter setup
    // [[codeblock]]
    // initialize (e.g. multi-threading backend)
    initialize(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);
    // [[/codeblock]]

    // Define the sub problem type tags (collection of model properties)
    using Model3D = Properties::TTag::TissueTransportModel;
    using Model1D = Properties::TTag::NetworkTransportModel;

    // Start the timer (to measure code execution time)
    Dune::Timer timer;

    // ### Grid setup
    // [[codeblock]]
    // Create grids (see input file groups "[Tissue.Grid]" and "[Network.Grid]")
    // See properties.hh for the types of Grid
    using Grid3D = GetPropType<Model3D, Properties::Grid>;
    GridManager<Grid3D> gridManager3D;
    gridManager3D.init("Tissue");

    using Grid1D = GetPropType<Model1D, Properties::Grid>;
    GridManager<Grid1D> gridManager1D;
    gridManager1D.init("Network");

    // get data from grid file like vessel radii and fluxes in the network
    auto gridData1D = gridManager1D.getGridData();

    // we compute on the leaf grid view (the refined elements)
    const auto& gridView3D = gridManager3D.grid().leafGridView();
    const auto& gridView1D = gridManager1D.grid().leafGridView();

    // create the finite volume grid geometries on top of the grid
    // this enables us to iterate over (sub-)control-volumes and faces
    using GridGeometry3D = GetPropType<Model3D, Properties::GridGeometry>;
    using GridGeometry1D = GetPropType<Model1D, Properties::GridGeometry>;
    auto gridGeometry3D = std::make_shared<GridGeometry3D>(gridView3D);
    auto gridGeometry1D = std::make_shared<GridGeometry1D>(gridView1D);

    // print info about the grid setup time
    const auto gridSetupTime = timer.elapsed();
    std::cout << "Setting up grid geometry took " << gridSetupTime << " seconds." << std::endl;
    // [[/codeblock]]
    //
    // ### Coupling manager
    //
    // [[codeblock]]
    // the mixed dimension type traits (create coupled model properties)
    using CoupledModelTraits = MultiDomainTraits<Model3D, Model1D>;
    constexpr auto domain3D = CoupledModelTraits::template SubDomain<0>::Index();
    constexpr auto domain1D = CoupledModelTraits::template SubDomain<1>::Index();

    // the coupling manager (handles data mapping between 1D and 3D domains)
    using CouplingManager = Properties::CouplingTransport;
    auto couplingManager = std::make_shared<CouplingManager>(gridGeometry3D, gridGeometry1D);

    // the problem (initial and boundary conditions, source and coupling terms)
    // these class are defined in `problem.hh` and then defines as "properties" in `properties.hh`
    // the constructor is whatever is defined as a constructor in `problem.hh`, so the user is
    // free to change this to whatever works. Here we added the argument `gridData1D` to hand
    // the grid parameters read from the grid file to the problem class instance
    using TransportProblem3D = GetPropType<Model3D, Properties::Problem>;
    using TransportProblem1D = GetPropType<Model1D, Properties::Problem>;
    auto transportProblem3D = std::make_shared<TransportProblem3D>(gridGeometry3D, couplingManager, "Tissue");
    auto transportProblem1D = std::make_shared<TransportProblem1D>(gridGeometry1D, couplingManager, gridData1D, "Network");

    // solve blood flow problem to obtain volume fluxes
    auto bloodVolumeFluxes = computeBloodVolumeFluxes(gridGeometry1D, gridData1D);
    transportProblem1D->spatialParams().setVolumeFluxes(bloodVolumeFluxes);
    // [[/codeblock]]
    //
    // ### Solution vector
    //
    // [[codeblock]]
    // the solution vector (here: mole fractions of the transported solute)
    CoupledModelTraits::SolutionVector sol;
    transportProblem3D->applyInitialSolution(sol[domain3D]);
    transportProblem1D->applyInitialSolution(sol[domain1D]);

    // this is a time-dependent problem, so we create a second solution vector that
    // represents the solution at the previous time level
    auto oldSol = sol;
    // [[/codeblock]]
    //
    // ### Coupling manager
    //
    // [[codeblock]]
    // initialize coupling
    couplingManager->init(transportProblem3D, transportProblem1D, sol);
    // in the case of 1D-3D the coupling needs point sources to be enabled for both sub-problems
    // each point sources constitutes an integration point for the coupling operator
    transportProblem3D->computePointSourceMap();
    transportProblem1D->computePointSourceMap();

    // the grid variables (secondary variables and caching of computed variables)
    // for instance the effective diffusion coefficient that is computed from
    // the free diffusion coefficient (fluid system) and the porosity and tortuosity (spatial parameter)
    using GridVariables3D = GetPropType<Model3D, Properties::GridVariables>;
    auto gridVariables3D = std::make_shared<GridVariables3D>(transportProblem3D, gridGeometry3D);
    gridVariables3D->init(sol[domain3D]);
    using GridVariables1D = GetPropType<Model1D, Properties::GridVariables>;
    auto gridVariables1D = std::make_shared<GridVariables1D>(transportProblem1D, gridGeometry1D);
    gridVariables1D->init(sol[domain1D]);
    // [[/codeblock]]
    //
    // ### VTK output
    //
    // [[codeblock]]
    // initialize the VTK output
    using Solution3D = std::decay_t<decltype(sol[domain3D])>;
    VtkOutputModule<GridVariables3D, Solution3D> vtkWriter3D(*gridVariables3D, sol[domain3D], transportProblem3D->name());
    GetPropType<Model3D, Properties::IOFields>::initOutputModule(vtkWriter3D);
    using Solution1D = std::decay_t<decltype(sol[domain1D])>;
    VtkOutputModule<GridVariables1D, Solution1D> vtkWriter1D(*gridVariables1D, sol[domain1D], transportProblem1D->name());
    GetPropType<Model1D, Properties::IOFields>::initOutputModule(vtkWriter1D);

    // we compute with mole fraction as primary variables because they are continuous across material interfaces
    // for the output, we compute fluid concentration c = rho*x and total/voxel concentration ct = porosity*c
    // mmol/l is equivalent to mol/m^3, the balance equation is in mol/s (change of tracer amount per time)
    vtkWriter3D.addVolumeVariable( [](const auto& v){ return v.molarDensity()*v.moleFraction(0,0); }, "c fluid (mmol/l)");
    vtkWriter3D.addVolumeVariable( [](const auto& v){ return v.molarDensity()*v.moleFraction(0,0)*v.porosity(); }, "c total (mmol/l)");

    // for the network domain we want to output the vessel lumen radius for visualization
    // and the outer vessel radius (where we couple with the tissue domain)
    vtkWriter1D.addField(transportProblem1D->spatialParams().getVesselRadii(), "vessel radius (m)");
    vtkWriter1D.addField(transportProblem1D->spatialParams().getOuterRadii(), "outer radius (m)");
    vtkWriter1D.addField(transportProblem1D->spatialParams().getVesselSegmentLengths(), "length");
    vtkWriter1D.addVolumeVariable( [](const auto& v){ return v.molarDensity()*v.moleFraction(0,0); }, "c fluid (mmol/l)");

    // write out initial solution
    vtkWriter3D.write(0.0, Dune::VTK::OutputType::appendedraw);
    vtkWriter1D.write(0.0, Dune::VTK::OutputType::appendedraw);
    // [[/codeblock]]
    //
    // ### Time loop configuration, assembler and linear solver
    //
    // [[codeblock]]
    ////////////////////////////////////////////////////////////
    // Time loop with assembler and linear solver
    ////////////////////////////////////////////////////////////

    // get some time loop parameters and make time loop
    const auto dt = getParam<double>("TimeLoop.DtInitial");
    const auto tEnd = getParam<double>("TimeLoop.TEnd");
    auto timeLoop = std::make_shared<CheckPointTimeLoop<double>>(0.0, dt, tEnd);

    // the assembler with time loop for a coupled transient problem
    using Assembler = MultiDomainFVAssembler<CoupledModelTraits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(transportProblem3D, transportProblem1D),
        std::make_tuple(gridGeometry3D, gridGeometry1D),
        std::make_tuple(gridVariables3D, gridVariables1D),
        couplingManager, timeLoop, oldSol
    );

    // construct the linear solver (see solver.hh)
    using LinearSolver = Example::BlockDiagILU0BiCGSTABSolver<typename Assembler::JacobianMatrix, typename Assembler::ResidualType>;
    auto linearSolver = std::make_shared<LinearSolver>();
    // [[/codeblock]]
    //
    // ### Output of the tracer concentration
    //
    // [[codeblock]]
    // output time and tracer amount in tissue
    const auto filePrefix = getParam<std::string>("Problem.Name");
    std::vector<Dune::FieldVector<double, 2>> tracerAmounts;
    tracerAmounts.reserve(std::size_t(std::ceil(tEnd/dt)));
    const auto initTracerAmount = transportProblem3D->computeTracerAmount(sol[domain3D], *gridVariables3D);
    std::cout << Fmt::format("Initial amount of tracer in tissue: {:.2e}\n", initTracerAmount);
    tracerAmounts.push_back(Dune::FieldVector<double, 2>{ 0.0, initTracerAmount });
    // [[/codeblock]]
    //
    // ### Time loop
    //
    // [[codeblock]]
    // start the time loop
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // construct some timers to measure the execution time of the time integration steps
        Dune::Timer assembleTimer(false), solveTimer(false), updateTimer(false), outputTimer(false);
        std::cout << "\nAssemble linear system " << std::endl;

        // assemble stiffness matrix and/or residual
        // we only assemble the Jacobian (stiffness matrix for a linear problem) once
        // in the subsequent time steps only the right-hand side (residual) changes
        assembleTimer.start();
        couplingManager->updateSolution(sol);
        if (timeLoop->timeStepIndex() == 0)
            assembler->assembleJacobianAndResidual(sol);
        else
            assembler->assembleResidual(sol);
        assembleTimer.stop();

        std::cout << "Solve linear system (" << linearSolver->name() << ")" << std::endl;

        // solve linear system
        // the first time we tell the solver about the Jacobian
        // then we solve with updated right-hand sides (residual)
        solveTimer.start();
        auto deltaSol = sol;
        if (timeLoop->timeStepIndex() == 0)
            linearSolver->setup(assembler->jacobian());

        const bool converged = linearSolver->solve(deltaSol, assembler->residual());
        if (!converged)
            DUNE_THROW(Dune::MathError, "Linear solver did not converge!");
        solveTimer.stop();

        // update solution and primary and secondary variables
        // and tell the coupling manager about the new solution
        updateTimer.start();
        sol -= deltaSol;
        couplingManager->updateSolution(sol);
        gridVariables3D->update(sol[domain3D]);
        gridVariables1D->update(sol[domain1D]);

        // make the new solution the old solution
        oldSol = sol;
        gridVariables3D->advanceTimeStep();
        gridVariables1D->advanceTimeStep();
        updateTimer.stop();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        outputTimer.start();
        vtkWriter3D.write(timeLoop->time(), Dune::VTK::OutputType::appendedraw);
        vtkWriter1D.write(timeLoop->time(), Dune::VTK::OutputType::appendedraw);

        // Compute total tracer mass in domain and write to file
        const auto tracerAmount = transportProblem3D->computeTracerAmount(sol[domain3D], *gridVariables3D);
        std::cout << Fmt::format("Current amount of tracer in ECS: {:.2e}\n", tracerAmount);
        tracerAmounts.push_back(Dune::FieldVector<double, 2>{ timeLoop->time(), tracerAmount });
        Dumux::writeContainerToFile(tracerAmounts, Fmt::format("{}_tracer_amounts.dat", filePrefix));
        outputTimer.stop();

        // output timings
        const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed() + outputTimer.elapsed();
        std::cout << "Assemble/solve/update/output time: "
                  <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                  <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                  <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)/"
                  <<  outputTimer.elapsed() << "(" << 100*outputTimer.elapsed()/elapsedTot << "%)"
                  << "\n";

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // use same time step size in next time step
        // this function can be used to implement time step adaptivity
        timeLoop->setTimeStepSize(dt);
    }
    // [[/codeblock]]
    //
    // ### Finalize
    //
    // [[codeblock]]
    // end the time loop and produce report
    timeLoop->finalize(gridView3D.comm());

    // output (manual) timings
    std::cout << "Transport computation took " << timer.elapsed() - gridSetupTime << " seconds." << std::endl;
    std::cout << "Simulation took " << timer.elapsed() << " seconds." << std::endl;

    // Print used/unused parameters at the end (good for debugging)
    Parameters::print();

    return 0;
}
// [[/codeblock]]
// [[/content]]
