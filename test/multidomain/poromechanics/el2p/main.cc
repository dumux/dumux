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

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

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

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////
    using TwoPTypeTag = Properties::TTag::TwoPSub;
    using PoroMechTypeTag = Properties::TTag::PoroElasticSub;

    // we simply extract the grid creator from one of the type tags
    using GridManager = Dumux::GridManager<GetPropType<TwoPTypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run stationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometries
    using TwoPFVGridGeometry = GetPropType<TwoPTypeTag, Properties::GridGeometry>;
    using PoroMechFVGridGeometry = GetPropType<PoroMechTypeTag, Properties::GridGeometry>;
    auto twoPFvGridGeometry = std::make_shared<TwoPFVGridGeometry>(leafGridView);
    auto poroMechFvGridGeometry = std::make_shared<PoroMechFVGridGeometry>(leafGridView);

    // the coupling manager
    using Traits = MultiDomainTraits<TwoPTypeTag, PoroMechTypeTag>;
    using CouplingManager = PoroMechanicsCouplingManager<Traits>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using TwoPProblem = GetPropType<TwoPTypeTag, Properties::Problem>;
    using PoroMechProblem = GetPropType<PoroMechTypeTag, Properties::Problem>;
    auto twoPSpatialParams = std::make_shared<typename TwoPProblem::SpatialParams>(twoPFvGridGeometry, couplingManager);
    auto twoPProblem = std::make_shared<TwoPProblem>(twoPFvGridGeometry, twoPSpatialParams, "TwoP");
    auto poroMechSpatialParams = std::make_shared<typename PoroMechProblem::SpatialParams>(poroMechFvGridGeometry, couplingManager);
    auto poroMechProblem = std::make_shared<PoroMechProblem>(poroMechFvGridGeometry, poroMechSpatialParams, "PoroElastic");

    // the solution vectors
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;

    static const auto twoPId = Traits::template SubDomain<0>::Index();
    static const auto poroMechId = Traits::template SubDomain<1>::Index();
    x[twoPId].resize(twoPFvGridGeometry->numDofs());
    x[poroMechId].resize(poroMechFvGridGeometry->numDofs());
    twoPProblem->applyInitialSolution(x[twoPId]);
    poroMechProblem->applyInitialSolution(x[poroMechId]);
    SolutionVector xOld = x;

    // initialize the coupling manager
    couplingManager->init(twoPProblem, poroMechProblem, x);

    // the grid variables
    using TwoPGridVariables = GetPropType<TwoPTypeTag, Properties::GridVariables>;
    using PoroMechGridVariables = GetPropType<PoroMechTypeTag, Properties::GridVariables>;
    auto twoPGridVariables = std::make_shared<TwoPGridVariables>(twoPProblem, twoPFvGridGeometry);
    auto poroMechGridVariables = std::make_shared<PoroMechGridVariables>(poroMechProblem, poroMechFvGridGeometry);
    twoPGridVariables->init(x[twoPId]);
    poroMechGridVariables->init(x[poroMechId]);

    // get some time loop parameters
    using Scalar = GetPropType<TwoPTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDT = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.Dt");

    // initialize the vtk output module
    using TwoPVtkOutputModule = Dumux::VtkOutputModule<TwoPGridVariables, GetPropType<TwoPTypeTag, Properties::SolutionVector>>;
    using PoroMechVtkOutputModule = Dumux::VtkOutputModule<PoroMechGridVariables, GetPropType<PoroMechTypeTag, Properties::SolutionVector>>;
    TwoPVtkOutputModule twoPVtkWriter(*twoPGridVariables, x[twoPId], twoPProblem->name());
    PoroMechVtkOutputModule poroMechVtkWriter(*poroMechGridVariables, x[poroMechId], poroMechProblem->name());

    // add output fields to writers
    using TwoPOutputFields = GetPropType<TwoPTypeTag, Properties::IOFields>;
    using PoroMechOutputFields = GetPropType<PoroMechTypeTag, Properties::IOFields>;
    TwoPOutputFields::initOutputModule(twoPVtkWriter);
    PoroMechOutputFields::initOutputModule(poroMechVtkWriter);

    // write initial solution
    twoPVtkWriter.write(0.0);
    poroMechVtkWriter.write(0.0);

    //instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDT);

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

    // the assembler
    using Assembler = Experimental::MultiStageMultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(twoPProblem, poroMechProblem),
                                                 std::make_tuple(twoPFvGridGeometry, poroMechFvGridGeometry),
                                                 std::make_tuple(twoPGridVariables, poroMechGridVariables),
                                                 couplingManager, timeSteppingMethod, xOld);

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

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();
        twoPGridVariables->advanceTimeStep();
        poroMechGridVariables->advanceTimeStep();

        // write vtk output
        twoPVtkWriter.write(timeLoop->time());
        poroMechVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        using TwoPPrimaryVariables = GetPropType<TwoPTypeTag, Properties::PrimaryVariables>;
        TwoPPrimaryVariables storage(0);
        const auto& twoPLocalResidual = assembler->localResidual(twoPId);
        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            const auto fvGeometry = localView(*twoPFvGridGeometry).bindElement(element);
            const auto elemVolVars = localView(twoPGridVariables->curGridVolVars()).bindElement(element, fvGeometry, x[twoPId]);
            storage += twoPLocalResidual.evalStorage(fvGeometry, elemVolVars)[0];
        }
        std::cout << "time, mass CO2 (kg), mass brine (kg):" << std::endl;
        std::cout << timeLoop->time() << " , " << storage[1] << " , " << storage[0] << std::endl;
        std::cout << "***************************************" << std::endl;

    } while (!timeLoop->finished());


    // output some Newton statistics
    nonLinearSolver->report();

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;

}
