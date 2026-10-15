// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/assembly/multistageassembler.hh>
#include <dumux/experimental/timestepping/generalizedalpha.hh>

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

    // get the refinement
    const auto refinement = getParam<int>("Grid.Refinement", 0);
    gridManager.grid().globalRefine(refinement);

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    const auto numElements = leafGridView.size(0);
    const auto numDofs = gridGeometry->numDofs();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    const std::string problemName = problem->name();
    const std::string csvFileName = problemName + "_dumux_box.csv";

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

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dtInitial = getParam<Scalar>("TimeLoop.DtInitial");
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dtInitial, tEnd);

    // the time integration method
    const auto method = [&]() -> std::shared_ptr<const Experimental::GeneralizedAlphaMethod<Scalar>>
    {
        using namespace Experimental::GeneralizedAlpha;
        const auto name = getParam<std::string>("TimeStepping.Method", "Newmark");
        if (name == "Newmark")
            return std::make_shared<Newmark<Scalar>>(getParam<Scalar>("TimeStepping.Beta", 0.25),
                                                     getParam<Scalar>("TimeStepping.Gamma", 0.5));

        const auto rhoInf = getParam<Scalar>("TimeStepping.SpectralRadius");
        if (name == "ChungHulbert")
            return std::make_shared<ChungHulbert<Scalar>>(rhoInf);
        else if (name == "HilberHughesTaylor")
            return std::make_shared<HilberHughesTaylor<Scalar>>(rhoInf);
        else if (name == "WoodBossakZienkiewicz")
            return std::make_shared<WoodBossakZienkiewicz<Scalar>>(rhoInf);
        DUNE_THROW(Dune::InvalidStateException, "Unknown time stepping method " << name);
    }();

    // the assembler for the implicit solves of the time integration method
    auto xOld = x;
    using Assembler = Experimental::MultiStageAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, method, xOld);

    // the linear solver
    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LAT>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver);

    // the beam starts at rest, and gravity acts from t = 0 on
    Experimental::GeneralizedAlphaTimeStepper<NewtonSolver> timeStepper(nonLinearSolver, method);
    if (getParam<bool>("TimeStepping.ConsistentInitialAcceleration", true))
    {
        // the undeformed beam is stress-free, so the free degrees of freedom accelerate with gravity
        SolutionVector velocity(x.size()), acceleration(x.size());
        velocity = 0.0;
        acceleration = 0.0;
        const auto gravity = getParam<Scalar>("Problem.Gravity");
        for (auto& a : acceleration)
            a[1] = -gravity;
        for (const auto& constraint : problem->constraints())
            acceleration[constraint.dofIndex()] = 0.0;
        timeStepper.initialize(x, 0.0, velocity, acceleration);
    }
    else
        timeStepper.initialize(x, 0.0);

    const int vtkInterval = getParam<int>("VTKOutput.Every", 5);

    // If the output already exist we overwrite it
    std::ofstream csvFile(csvFileName, std::ios::out | std::ios::trunc);
    csvFile << "t,dt,level,nel,ndof,ux,uy\n";

    // time loop
    timeLoop->start(); do
    {
        // time integration
        timeStepper.step(x, timeLoop->time(), timeLoop->timeStepSize());
        xOld = x;

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        const auto displacementA = problem->evalControlPointDisplacement(*gridVariables, x);
        csvFile << timeLoop->time() << "," << timeLoop->timeStepSize() << ","
                << refinement << "," <<numElements << "," << numDofs << ","
                << displacementA[0] << "," << displacementA[1] << "\n";

        // write VTK output
        if (timeLoop->timeStepIndex() % vtkInterval == 0 || timeLoop->finished())
            vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    return 0;
}
