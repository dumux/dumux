// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoromechanicsTests
 * \brief Test for a single-phase elastic coupled model.
 */

#include <config.h>
#include <iostream>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/experimental/assembly/multistagemultidomainfvassembler.hh>
#include <dumux/experimental/timestepping/timelevel.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>
#include <dumux/discretization/evalgradients.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameter tree
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////
    using OnePTypeTag = Properties::TTag::OnePSub;
    using PoroMechTypeTag = Properties::TTag::PoroElasticSub;

    // we simply extract the grid creator from one of the type tags
    using GridManager = Dumux::GridManager<GetPropType<OnePTypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run stationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometries
    using OnePFVGridGeometry = GetPropType<OnePTypeTag, Properties::GridGeometry>;
    using PoroMechFVGridGeometry = GetPropType<PoroMechTypeTag, Properties::GridGeometry>;
    auto onePFvGridGeometry = std::make_shared<OnePFVGridGeometry>(leafGridView);
    auto poroMechFvGridGeometry = std::make_shared<PoroMechFVGridGeometry>(leafGridView);

    // the coupling manager
    using Traits = MultiDomainTraits<OnePTypeTag, PoroMechTypeTag>;
    using CouplingManager = PoroMechanicsCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // create analytical solution
    using Scalar = GetPropType<OnePTypeTag, Properties::Scalar>;
    auto mandelAnalyticalSolution = std::make_shared<MandelAnalyticalSolution<Scalar>>();

    // print beta values to file
    mandelAnalyticalSolution->printBetaValues();

    // the problems (boundary conditions)
    using OnePProblem = GetPropType<OnePTypeTag, Properties::Problem>;
    using PoroMechProblem = GetPropType<PoroMechTypeTag, Properties::Problem>;
    auto onePSpatialParams = std::make_shared<typename OnePProblem::SpatialParams>(onePFvGridGeometry, couplingManager);
    auto onePProblem = std::make_shared<OnePProblem>(onePFvGridGeometry, onePSpatialParams, "OneP");
    auto poroMechSpatialParams = std::make_shared<typename PoroMechProblem::SpatialParams>(poroMechFvGridGeometry, couplingManager);
    auto poroMechProblem = std::make_shared<PoroMechProblem>(poroMechFvGridGeometry, poroMechSpatialParams, "PoroElastic");
    poroMechSpatialParams->setAnalyticalSolution(*mandelAnalyticalSolution);
    onePSpatialParams->setAnalyticalSolution(*mandelAnalyticalSolution);

    // the solution vectors
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;

    static const auto onePId = Traits::template SubDomain<0>::Index();
    static const auto poroMechId = Traits::template SubDomain<1>::Index();
    x[onePId].resize(onePFvGridGeometry->numDofs());
    x[poroMechId].resize(poroMechFvGridGeometry->numDofs());
    onePProblem->applyInitialSolution(x[onePId]);
    poroMechProblem->applyInitialSolution(x[poroMechId]);
    SolutionVector xOld = x;

    // vectors for the exact solution
    auto pExact = x[onePId];
    auto uExact = x[poroMechId];


    // initialize the coupling manager
    couplingManager->init(onePProblem, poroMechProblem, x);

    // the grid variables
    using OnePGridVariables = GetPropType<OnePTypeTag, Properties::GridVariables>;
    using PoroMechGridVariables = GetPropType<PoroMechTypeTag, Properties::GridVariables>;
    auto onePGridVariables = std::make_shared<OnePGridVariables>(onePProblem, onePFvGridGeometry);
    auto poroMechGridVariables = std::make_shared<PoroMechGridVariables>(poroMechProblem, poroMechFvGridGeometry);
    onePGridVariables->init(x[onePId]);
    poroMechGridVariables->init(x[poroMechId]);

    // get some time loop parameters
    using Scalar = GetPropType<OnePTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDT = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.Dt");

    // initialize the vtk output module
    using OnePVtkOutputModule = Dumux::VtkOutputModule<OnePGridVariables, GetPropType<OnePTypeTag, Properties::SolutionVector>>;
    using PoroMechVtkOutputModule = Dumux::VtkOutputModule<PoroMechGridVariables, GetPropType<PoroMechTypeTag, Properties::SolutionVector>>;
    OnePVtkOutputModule onePVtkWriter(*onePGridVariables, x[onePId], onePProblem->name());
    onePVtkWriter.addVolumeVariable([](const auto &v) { return v.porosity(); }, "porosity");
    onePVtkWriter.addVolumeVariable([](const auto &v) { return v.density(); }, "rho");
    onePVtkWriter.addVolumeVariable([](const auto &v) { return v.permeability(); }, "Permeability");

    PoroMechVtkOutputModule poroMechVtkWriter(*poroMechGridVariables, x[poroMechId], poroMechProblem->name());

    // add output fields to writers
    using OnePOutputFields = GetPropType<OnePTypeTag, Properties::IOFields>;
    using PoroMechOutputFields = GetPropType<PoroMechTypeTag, Properties::IOFields>;
    OnePOutputFields::initOutputModule(onePVtkWriter);
    PoroMechOutputFields::initOutputModule(poroMechVtkWriter);

    // compute divergence of displacement field for output as cell data
    std::vector<Scalar> divUField(onePFvGridGeometry->gridView().size(0));  // one value per cell
    for (const auto& element : elements(onePFvGridGeometry->gridView()))
    {
        const auto elemGeom = element.geometry();
        const auto eIdx = onePFvGridGeometry->elementMapper().index(element);

        // get corresponding poromechanics element solution
        auto poroMechElemSol = elementSolution(element, x[poroMechId], *poroMechFvGridGeometry);

        // compute divergence of displacement at element center
        const auto gradU = Dumux::evalGradients(element, elemGeom, *poroMechFvGridGeometry, poroMechElemSol, elemGeom.center());

        Scalar divU = 0.0;
        for (int dir = 0; dir < PoroMechGridVariables::GridGeometry::GridView::dimension; ++dir)
            divU += gradU[dir][dir];

        divUField[eIdx] = divU;  // store scalar value per cell
    }

    // add divergence of displacement output to one-phase writer as cell data
    onePVtkWriter.addField(divUField, "divU", Dumux::Vtk::FieldType::element);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDT);

    // get the corner positions for the analytical solution
    Dune::BlockVector<Dune::FieldVector<Scalar, 2>> cornerPos;
    cornerPos.resize(uExact.size());
    for(const auto& vertex: vertices(poroMechFvGridGeometry->gridView()))
    {
        const auto vIdx = poroMechFvGridGeometry->vertexMapper().index(vertex);
        cornerPos[vIdx] = vertex.geometry().center();
    }

    Dumux::parallelFor(onePFvGridGeometry->gridView().size(0), [&](const std::size_t eIdx)
    { pExact[eIdx] = mandelAnalyticalSolution->initialPressure(onePFvGridGeometry->element(eIdx).geometry().center()); });

    Dumux::parallelFor(uExact.size(), [&](const std::size_t vIdx)
    { uExact[vIdx] = mandelAnalyticalSolution->initialDisplacement(cornerPos[vIdx]); });

    onePVtkWriter.addField(pExact, "pExact");
    poroMechVtkWriter.addField(uExact, "uExact");

    // write initial solution
    onePVtkWriter.write(0.0);
    poroMechVtkWriter.write(0.0);

    const auto timeSteppingMethod = []() -> std::shared_ptr<Experimental::MultiStageMethod<double>>
    {
        const auto timeStepScheme = getParam<std::string>("TimeLoop.Scheme", "DIRK3");
        if (timeStepScheme == "DIRK3")
            return std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<double>>();
        else if (timeStepScheme == "ImplicitEuler")
            return std::make_shared<Experimental::MultiStage::ImplicitEuler<double>>();
        else if (timeStepScheme == "CrankNicolson")
            return std::make_shared<Experimental::MultiStage::Theta<double>>(0.5);
        else
            DUNE_THROW(Dumux::ParameterException,
                "Unknown time stepping scheme. Possible values are "
                << "DIRK3 or ImplicitEuler or CrankNicolson");
    }();

    poroMechProblem->setTimeLoop(timeLoop);

    // the assembler
    using Assembler = Experimental::MultiStageMultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(onePProblem, poroMechProblem),
        std::make_tuple(onePFvGridGeometry, poroMechFvGridGeometry),
        std::make_tuple(onePGridVariables, poroMechGridVariables),
        couplingManager, timeSteppingMethod, xOld
    );

    // the linear solver
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    using TimeStepper = Experimental::MultiStageTimeStepper<NewtonSolver>;
    TimeStepper timeStepper(nonLinearSolver, timeSteppingMethod);

    // time loop
    timeLoop->start(); do
    {
        // solve the non-linear system with time step control
        timeStepper.step(x, timeLoop->time(), timeLoop->timeStepSize());

        // make the new solution the old solution
        xOld = x;

        // advance the time loop to the next step
        timeLoop->advanceTimeStep();
        onePGridVariables->advanceTimeStep();
        poroMechGridVariables->advanceTimeStep();

        Dumux::parallelFor(onePFvGridGeometry->gridView().size(0), [&](const std::size_t eIdx)
        { pExact[eIdx] = mandelAnalyticalSolution->pressure(onePFvGridGeometry->element(eIdx).geometry().center(), timeLoop->time()); });

        Dumux::parallelFor(uExact.size(), [&](const std::size_t vIdx)
        { uExact[vIdx] = mandelAnalyticalSolution->displacement(cornerPos[vIdx], timeLoop->time()); });

        // update divergence of displacement field for current time step
        for (const auto& element : elements(onePFvGridGeometry->gridView()))
        {
            const auto elemGeom = element.geometry();
            const auto eIdx = onePFvGridGeometry->elementMapper().index(element);

            // get corresponding poromechanics element solution
            auto poroMechElemSol = elementSolution(element, x[poroMechId], *poroMechFvGridGeometry);

            // compute divergence of displacement at element center
            const auto gradU = Dumux::evalGradients(element, elemGeom, *poroMechFvGridGeometry, poroMechElemSol, elemGeom.center());

            Scalar divU = 0.0;
            for (int dir = 0; dir < PoroMechGridVariables::GridGeometry::GridView::dimension; ++dir)
                divU += gradU[dir][dir];

            divUField[eIdx] = divU;  // update scalar value per cell
        }

        // write vtk output
        onePVtkWriter.write(timeLoop->time());
        poroMechVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

    } while (!timeLoop->finished());

    // output some Newton statistics
    nonLinearSolver->report();

    timeLoop->finalize(leafGridView.comm());

    if (leafGridView.comm().rank() == 0)
        Parameters::print();

    return 0;
}
