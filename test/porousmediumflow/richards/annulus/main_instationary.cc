// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/integrate.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/io/format.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // We parse the command line arguments.
    Parameters::init(argc, argv);

    // Initialize timer
    Dune::Timer timer;

    // Convenience alias for the type tag of the problem.
    using TypeTag = Properties::TTag::RichardsAnnulus;

    // The grid manager can be used to create a grid from the input file
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    // We compute on the leaf grid view.
    const auto& leafGridView = gridManager.grid().leafGridView();

    // instantiate the grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // Initialize the problem and grid variables
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);
    problem->disableSource();
    problem->disableStationaryInitialSolution();
    problem->enableInnerNeumannBC();

    // We define a function to update the discrete analytical solution vector
    // using the exactSolution() function in the problem
    // Note that this refers to the exact solution in the stationary case
    // In the instationary case, this analytical solution is just an approximation
    const auto updateAnalyticalSolution = [&](auto& pAnalyticApproximation, const double innerPsi)
    {
        pAnalyticApproximation.resize(gridGeometry->numDofs());
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
                pAnalyticApproximation[scv.dofIndex()] = problem->exactSolution(scv.dofPosition(), innerPsi);
        }
    };

    // instantiate and initialize the discrete and exact solution vectors
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector p(gridGeometry->numDofs()); problem->applyInitialSolution(p);
    const auto initialPressure = problem->headToPressure(getParam<double>("Problem.InitialPressureHeadInCm", -100));
    SolutionVector pAnalyticApproximation; updateAnalyticalSolution(pAnalyticApproximation, problem->pressureToPsi(initialPressure));
    auto pOld = p;

    // instantiate and initialize the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(p);

    auto gridVariablesAnalytic = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariablesAnalytic->init(pAnalyticApproximation);
    const auto initialSaturation = problem->saturation(initialPressure);

    // get some time loop parameters
    const auto tEnd = getParam<double>("TimeLoop.TEnd", 24*60*60*40);
    const auto maxDt = getParam<double>("TimeLoop.MaxTimeStepSize", 0.5*24*60*60);
    auto dt = getParam<double>("TimeLoop.DtInitial", 0.5*24*60*60);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<double>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    std::vector<double> storageDerivative(leafGridView.size(Grid::dimension), 0.0);

    // initialize VTK output
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, p, problem->name());
    GetPropType<TypeTag, Properties::IOFields>::initOutputModule(vtkWriter);
    vtkWriter.addField(pAnalyticApproximation, "pAnalyticApproximation");
    vtkWriter.addField(storageDerivative, "dS/dt");
    vtkWriter.addVolumeVariable([&](const auto& v){
        return (v.pressure(0) - problem->nonwettingReferencePressure())*100/(v.density(0)*9.81); },
        "head in cm"
    );
    vtkWriter.write(0.0);

    // instantiate the solver
    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, pOld);

    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    // time loop
    bool finished = false;
    const auto wiltingPressure = problem->headToPressure(-15000);
    timeLoop->start(); do
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(p, *timeLoop);

        auto interfacePressure = p[0][0];
        while (interfacePressure < 1.001*wiltingPressure || interfacePressure > 1.1e5)
        {
            p = pOld;
            assembler->resetTimeStep(p);
            const auto dt = timeLoop->timeStepSize() * 0.5;
            std::cout << Fmt::format("Retry with time step size: {} s", dt) << std::endl;
            timeLoop->setTimeStepSize(dt);
            nonLinearSolver.solve(p, *timeLoop);
            interfacePressure = p[0][0];
        }

        std::cout << Fmt::format("Root-soil interface pressure: {} Pa", interfacePressure) << std::endl;
        problem->updateStorageDerivative(p, pOld, *gridVariables, timeLoop->timeStepSize(), storageDerivative);

        // make the new solution the old solution
        pOld = p;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // compute and compare with analytic approximation
        updateAnalyticalSolution(pAnalyticApproximation, problem->pressureToPsi(interfacePressure));
        gridVariablesAnalytic->update(pAnalyticApproximation);
        const auto timePrediction = problem->computeTime(pAnalyticApproximation, *gridVariablesAnalytic, initialSaturation);
        std::cout << Fmt::format("It took {} days until this pressure was reached\n", timeLoop->time()/86400.0);
        std::cout << Fmt::format("The analytical approximation predicted {} days\n", timePrediction/86400.0);

        // write vtk output
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        // check if we are done
        if (interfacePressure <= wiltingPressure)
            finished = true;

        finished |= timeLoop->finished();

    } while (!finished);

    const auto timePrediction = problem->computeTime(pAnalyticApproximation, *gridVariablesAnalytic, initialSaturation);
    if (auto d = std::abs(timeLoop->time()-timePrediction)/timePrediction; d > 0.05)
        DUNE_THROW(Dune::Exception, "Deviation between estimate and numerical solution too big: " << d*100 << "%");

    timeLoop->finalize(leafGridView.comm());
    nonLinearSolver.report();
    Parameters::print();

    return 0;
}
