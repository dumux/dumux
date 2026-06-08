// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>
#include <fstream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/geometry/diameter.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // define the sub problem type tags
    using BulkTypeTag = Properties::TTag::Soil;
    using LowDimTypeTag = Properties::TTag::Voids;

    // initialize the grid for the 3D domain
    using BulkGridManager = Dumux::GridManager<GetPropType<BulkTypeTag, Properties::Grid>>;
    BulkGridManager bulkGridManager;
    bulkGridManager.init("Soil"); // pass parameter group

    // initialize the grid for the 1D domain
    using LowDimGridManager = Dumux::GridManager<GetPropType<LowDimTypeTag, Properties::Grid>>;
    LowDimGridManager lowDimGridManager;
    lowDimGridManager.init("Voids"); // pass parameter group

    // create 3D grid geometry and grid view
    using BulkFVGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    const auto& bulkGridView = bulkGridManager.grid().leafGridView();
    auto bulkFvGridGeometry = std::make_shared<BulkFVGridGeometry>(bulkGridView);

    // create 1D grid geometry and grid view
    using LowDimFVGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
    const auto& lowDimGridView = lowDimGridManager.grid().leafGridView();
    auto lowDimFvGridGeometry = std::make_shared<LowDimFVGridGeometry>(lowDimGridView);

    // the mixed dimension type traits
    using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    constexpr auto bulkIdx = Traits::template SubDomain<0>::Index();
    constexpr auto lowDimIdx = Traits::template SubDomain<1>::Index();

    // the coupling manager
    using CouplingManager = GetPropType<BulkTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(bulkFvGridGeometry, lowDimFvGridGeometry);
    auto bulkExchangeFluxCalculator = std::make_shared<ExchangeFluxCalculator<BulkTypeTag>>(couplingManager);
    auto lowDimExchangeFluxCalculator = std::make_shared<ExchangeFluxCalculator<LowDimTypeTag>>(couplingManager);

    // the low dim spatial parameters
    using LowDimSpatialParams = GetPropType<LowDimTypeTag, Properties::SpatialParams>;
    auto lowDimSpatialParams = std::make_shared<LowDimSpatialParams>(lowDimFvGridGeometry, lowDimGridManager.getGridData());
    // the bulk spatial parameters
    using BulkSpatialParams = GetPropType<BulkTypeTag, Properties::SpatialParams>;
    auto bulkSpatialParams = std::make_shared<BulkSpatialParams>(bulkFvGridGeometry);

    // the problem (initial and boundary conditions)
    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, bulkSpatialParams, couplingManager, bulkExchangeFluxCalculator);
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, lowDimSpatialParams, couplingManager, lowDimExchangeFluxCalculator);

    // the solution vector
    Traits::SolutionVector sol;
    sol[bulkIdx].resize(bulkFvGridGeometry->numDofs());
    sol[lowDimIdx].resize(lowDimFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(sol[bulkIdx]);
    lowDimProblem->applyInitialSolution(sol[lowDimIdx]);
    auto oldSol = sol;

    couplingManager->init(bulkProblem, lowDimProblem, sol);
    bulkProblem->computePointSourceMap();
    lowDimProblem->computePointSourceMap();

    // the grid variables
    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    bulkGridVariables->init(sol[bulkIdx]);
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    lowDimGridVariables->init(sol[lowDimIdx]);

    // // pass gridVars to coupling manager
    // couplingManager->setGridVariables(std::make_tuple(bulkGridVariables, lowDimGridVariables));

    // get some time loop parameters
    using Scalar = Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto episodeLength = getParam<Scalar>("TimeLoop.EpisodeLength");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize", episodeLength);
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // initialize the vtk output module
    using BulkSolutionVector = std::decay_t<decltype(sol[bulkIdx])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, sol[bulkIdx], bulkProblem->name());
    GetPropType<BulkTypeTag, Properties::IOFields>::initOutputModule(bulkVtkWriter);
    bulkProblem->addVtkOutputFields(bulkVtkWriter);
    bulkVtkWriter.write(0.0, Dune::VTK::base64);

    using LowDimSolutionVector = std::decay_t<decltype(sol[lowDimIdx])>;
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, sol[lowDimIdx], lowDimProblem->name());
    GetPropType<LowDimTypeTag, Properties::IOFields>::initOutputModule(lowDimVtkWriter);
    using LowDimVelocityOutput = GetPropType<LowDimTypeTag, Properties::VelocityOutput>;
    lowDimVtkWriter.addVelocityOutput(std::make_shared<LowDimVelocityOutput>(*lowDimGridVariables));
    lowDimProblem->addVtkOutputFields(lowDimVtkWriter);
    lowDimVtkWriter.write(0.0, Dune::VTK::base64);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                 std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                 std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                 couplingManager, timeLoop, oldSol);

    // the linear solver
    // using LinearSolver = BlockDiagILU0BiCGSTABSolver;
    using LinearSolver = BlockDiagAMGBiCGSTABSolver;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // initialize the custom output to the csv file
    lowDimProblem->initializeOutput();

    // time loop
    timeLoop->setPeriodicCheckPoint(episodeLength);
    timeLoop->start(); do
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(sol, *timeLoop);

        // make the new solution the old solution
        oldSol = sol;
        bulkGridVariables->advanceTimeStep();
        lowDimGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        if (timeLoop->isCheckPoint())
        {
            // output the source terms
            bulkProblem->computeSourceIntegral(sol[bulkIdx], *bulkGridVariables);
            lowDimProblem->computeSourceIntegral(sol[lowDimIdx], *lowDimGridVariables);

            lowDimProblem->updateNusseltNumbers();

            bulkVtkWriter.write(timeLoop->time(), Dune::VTK::base64);
            lowDimVtkWriter.write(timeLoop->time(), Dune::VTK::base64);

            lowDimProblem->writeOutput(lowDimGridVariables,timeLoop, sol[lowDimIdx]);
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(mpiHelper.getCommunication());

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
} // end main
