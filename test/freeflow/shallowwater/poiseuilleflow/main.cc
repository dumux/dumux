// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_ug.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // the property tag
    using TypeTag = Properties::TTag::PoiseuilleFlow;

    // create grid
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // problem, in which we define the boundary and initial conditions.
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // We then initialize the predefined model-specific output vtk output.
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables,x, problem->name());
    vtkWriter.addField(problem->getExactWaterDepth(), "exactWaterDepth");
    vtkWriter.addField(problem->getExactVelocityX(), "exactVelocityX");
    vtkWriter.addField(problem->getExactVelocityY(), "exactVelocityY");
    problem->updateAnalyticalSolution();
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    // instantiate time loop
    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    if (hasParam("TimeLoop.PrintoutTimes"))
    {
        const auto timeLoopPrintoutTimes = getParam<std::vector<double>>("TimeLoop.PrintoutTimes");
        timeLoop->setCheckPoint(timeLoopPrintoutTimes.begin(), timeLoopPrintoutTimes.end());
    }

    timeLoop->setCheckPoint(tEnd);

    // assembler & solver
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);

        // update the analytical solution
        problem->updateAnalyticalSolution();

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        if (timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));


    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
