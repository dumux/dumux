// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include "analyticalsolution.hh"
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
    using BulkTypeTag = Properties::TTag::Rock;
    using LowDimTypeTag = Properties::TTag::Well;

    // initialize the grid for the 3D domain
    using BulkGridManager = Dumux::GridManager<GetPropType<BulkTypeTag, Properties::Grid>>;
    BulkGridManager bulkGridManager;
    bulkGridManager.init("Rock"); // pass parameter group

    // initialize the grid for the 1D domain
    using LowDimGridManager = Dumux::GridManager<GetPropType<LowDimTypeTag, Properties::Grid>>;
    LowDimGridManager lowDimGridManager;
    lowDimGridManager.init("Well"); // pass parameter group

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

    // get some time loop parameters
    using Scalar = Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto episodeLength = getParam<Scalar>("TimeLoop.EpisodeLength");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize", episodeLength);
    // a smaller maximum time step size until the time EarlyPeriod, which resolves the flushing of
    // the initial pipe content and the fast early cooling of the rock around the wellbore
    const auto earlyPeriod = getParam<Scalar>("TimeLoop.EarlyPeriod", 0.0);
    const auto maxDtEarly = getParam<Scalar>("TimeLoop.MaxTimeStepSizeEarly", maxDt);
    const auto maxTimeStepSize = [&](const Scalar time)
    { return time < earlyPeriod ? maxDtEarly : maxDt; };
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // initialize the vtk output module
    using BulkSolutionVector = std::decay_t<decltype(sol[bulkIdx])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, sol[bulkIdx], bulkProblem->name());
    GetPropType<BulkTypeTag, Properties::IOFields>::initOutputModule(bulkVtkWriter);
    bulkProblem->addVtkOutputFields(bulkVtkWriter);
    // the bulk output is large on fine grids, the convergence test switches it off
    const bool enableBulkOutput = getParam<bool>("Vtk.EnableBulkOutput", true);
    if (enableBulkOutput)
        bulkVtkWriter.write(0.0);

    using LowDimSolutionVector = std::decay_t<decltype(sol[lowDimIdx])>;
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, sol[lowDimIdx], lowDimProblem->name());
    GetPropType<LowDimTypeTag, Properties::IOFields>::initOutputModule(lowDimVtkWriter);
    using LowDimVelocityOutput = GetPropType<LowDimTypeTag, Properties::VelocityOutput>;
    lowDimVtkWriter.addVelocityOutput(std::make_shared<LowDimVelocityOutput>(*lowDimGridVariables));
    lowDimProblem->addVtkOutputFields(lowDimVtkWriter);
    lowDimVtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxTimeStepSize(0.0));

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                 std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                 std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                 couplingManager, timeLoop, oldSol);

    // the linear solver
    using LinearSolver = BlockDiagAMGBiCGSTABSolver;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // the reference solution, the radii are the same in all elements of the well
    using Reference = WellboreReferenceSolution<Scalar>;
    using Variant = Reference::Variant;
    const Reference reference(lowDimSpatialParams->innerRadius(0), lowDimSpatialParams->radius(0));
    const auto referenceCelsius = [&](const Scalar l, const Scalar t, const Variant variant)
    { return reference.temperature(l, t, variant) - 273.15; };

    // the distance from the inlet at the top and the index of all vertices of the well, sorted
    // along the well, so that the last one is the outlet
    std::vector<std::pair<Scalar, std::size_t>> wellDofs(lowDimFvGridGeometry->numDofs());
    for (const auto& element : elements(lowDimGridView))
    {
        const auto fvGeometry = localView(*lowDimFvGridGeometry).bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
            wellDofs[scv.dofIndex()] = {lowDimFvGridGeometry->bBoxMax()[2] - scv.dofPosition()[2], scv.dofIndex()};
    }
    std::sort(wellDofs.begin(), wellDofs.end());
    const Scalar outletDistance = wellDofs.back().first;
    const std::size_t outletDof = wellDofs.back().second;

    using LowDimIndices = GetPropType<LowDimTypeTag, Properties::ModelTraits>::Indices;
    const auto wellCelsius = [&](const std::size_t dofIdx)
    { return sol[lowDimIdx][dofIdx][LowDimIndices::temperatureIdx] - 273.15; };

    // the reference solutions at the outlet, independent of the output interval and dense
    // around the arrival of the first inlet water, where the outlet temperature has a kink
    const std::string outputName = getParam<std::string>("Vtk.OutputName");
    std::ofstream outletReferenceFile(outputName + "_outlet_reference.csv");
    outletReferenceFile << "time [d], T reference (a) [°C], T reference (b) [°C]\n" << std::setprecision(15);
    const auto writeOutletReference = [&](const Scalar t)
    {
        outletReferenceFile << t/86400 << "," << referenceCelsius(outletDistance, t, Variant::benchmark)
                            << "," << referenceCelsius(outletDistance, t, Variant::embedded) << "\n";
    };
    const Scalar denseEnd = 4.0*reference.travelTime(outletDistance);
    for (int i = 1; i <= 200; ++i)
        writeOutletReference(i*denseEnd/200);
    for (int i = 1; i <= 400; ++i)
        writeOutletReference(denseEnd + i*(tEnd - denseEnd)/400);
    outletReferenceFile.close();

    // the outlet temperature over time. The errors are measured against the embedded variant (b),
    // the exact solution of the continuous problem that the surface coupling discretizes
    std::ofstream outletFile(lowDimProblem->name() + ".csv");
    outletFile << "time [d], T [°C], T reference (a) [°C], T reference (b) [°C]\n" << std::setprecision(15);
    Scalar maxOutletError = 0.0;
    Scalar maxOutletErrorAfterOneDay = 0.0;

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

            if (enableBulkOutput)
                bulkVtkWriter.write(timeLoop->time());
            lowDimVtkWriter.write(timeLoop->time());

            // compare the outlet temperature with the reference solution
            const Scalar time = timeLoop->time();
            const Scalar outletTemperature = wellCelsius(outletDof);
            const Scalar outletReference = referenceCelsius(outletDistance, time, Variant::embedded);
            const Scalar outletError = std::abs(outletTemperature - outletReference);
            maxOutletError = std::max(maxOutletError, outletError);
            if (time >= 86400.0)
                maxOutletErrorAfterOneDay = std::max(maxOutletErrorAfterOneDay, outletError);

            outletFile << time/86400 << "," << outletTemperature << ","
                       << referenceCelsius(outletDistance, time, Variant::benchmark) << ","
                       << outletReference << std::endl;
            std::cout << std::setprecision(8) << "Outlet temperature [°C]: " << outletTemperature
                      << ", reference (b): " << outletReference << '\n';
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setMaxTimeStepSize(maxTimeStepSize(timeLoop->time()));
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(mpiHelper.getCommunication());

    // the reference solutions along the well at the final time, independent of the grid
    std::ofstream profileReferenceFile(outputName + "_profile_reference.csv");
    profileReferenceFile << "distance from the inlet [m], T reference (a) [°C], T reference (b) [°C]\n" << std::setprecision(15);
    for (int i = 0; i <= 200; ++i)
    {
        const Scalar distance = i*outletDistance/200;
        profileReferenceFile << distance << "," << referenceCelsius(distance, tEnd, Variant::benchmark)
                             << "," << referenceCelsius(distance, tEnd, Variant::embedded) << "\n";
    }
    profileReferenceFile.close();

    // compare the temperature along the well at the final time with the reference solution
    std::ofstream profileFile(outputName + "_profile.csv");
    profileFile << "distance from the inlet [m], T [°C], T reference (a) [°C], T reference (b) [°C]\n" << std::setprecision(15);
    Scalar maxProfileError = 0.0;
    for (const auto& [distance, dofIdx] : wellDofs)
    {
        const Scalar temperature = wellCelsius(dofIdx);
        const Scalar profileReference = referenceCelsius(distance, tEnd, Variant::embedded);
        maxProfileError = std::max(maxProfileError, std::abs(temperature - profileReference));
        profileFile << distance << "," << temperature << ","
                    << referenceCelsius(distance, tEnd, Variant::benchmark) << "," << profileReference << "\n";
    }

    std::cout << std::setprecision(8)
              << "\nMaximum error with regard to the reference solution (b) [K]\n"
              << "  outlet temperature, all outputs:        " << maxOutletError << "\n"
              << "  outlet temperature, outputs from 1 d on: " << maxOutletErrorAfterOneDay << "\n"
              << "  temperature along the well at the end:  " << maxProfileError << "\n";

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
