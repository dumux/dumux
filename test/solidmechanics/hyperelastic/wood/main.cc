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
#include <dumux/common/timeloop.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::WoodTest;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(getParam<Scalar>("TimeLoop.MaxTimeStepSize"));

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    vtkWriter.addVolumeVariable([](const auto& v){
        return Dune::FieldVector<double, 2>{v.displacement(0), v.displacement(1)};
    }, "d");
    vtkWriter.addVolumeVariable([](const auto& v){ return v.moistureContent(); }, "m");

    // local material frame (constant in time): radial and tangential basis vectors per element
    std::vector<Dune::FieldVector<double, 2>> eR(gridGeometry->gridView().size(0));
    std::vector<Dune::FieldVector<double, 2>> eT(gridGeometry->gridView().size(0));
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        const auto eIdx = gridGeometry->elementMapper().index(element);
        const auto Q = problem->spatialParams().rotation(element.geometry().center());
        eR[eIdx] = {Q[0][0], Q[1][0]}; // first column of Q
        eT[eIdx] = {Q[0][1], Q[1][1]}; // second column of Q
    }
    vtkWriter.addField(eR, "eR");
    vtkWriter.addField(eT, "eT");

    vtkWriter.write(0.0);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LAT>;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver);

    const int vtkInterval = getParam<int>("VTKOutput.Every", 1);

    timeLoop->start(); do
    {
        problem->setTime(timeLoop->time() + timeLoop->timeStepSize());
        nonLinearSolver->solve(x, *timeLoop);

        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        if (timeLoop->timeStepIndex() % vtkInterval == 0 || timeLoop->finished())
            vtkWriter.write(timeLoop->time());

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver->suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    return 0;
}
