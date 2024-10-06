// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <iostream>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/experimental/timestepping/newmarkbeta.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::DynamicHyperelasticityTest;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameter tree
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // initialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    vtkWriter.addVolumeVariable([](const auto& v){
        return Dune::FieldVector<double, 2>{v.displacement(0), v.displacement(1)};
    }, "d");
    vtkWriter.write(0.0);

    // solution at previous time
    auto xOld = x;

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);

    // the time stepping scheme
    auto newmarkBeta = std::make_shared<Experimental::NewmarkBeta<Scalar, SolutionVector>>();
    newmarkBeta->initialize(x);
    problem->setNewmarkScheme(newmarkBeta);

    // the assembler with time loop for a transient problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LAT>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver);

    const int vtkInterval = getParam<int>("VTKOutput.Every", 5);

    // time loop
    timeLoop->start(); do
    {
        // linearize & solve
        nonLinearSolver->solve(x);

        // update the solution in the time stepping scheme
        newmarkBeta->update(dt, x);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write VTK output
        if (timeLoop->timeStepIndex() % vtkInterval == 0 || timeLoop->finished())
            vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    return 0;
}
