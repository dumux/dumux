// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \author Timo Koch <timokoch@uio.no>
 * \brief Benchmark cases C1.2a/b
 * Cases are described in Schnepf et al (2020) https://doi.org/10.3389/fpls.2020.00316
 * We use the coupling method described in Koch et al (2022), https://doi.org/10.1016/j.jcp.2021.110823
 */
#include <config.h>

#include <iostream>
#include <fstream>
#include <algorithm>

#include <dune/common/timer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/initialize.hh>
#include <dumux/geometry/diameter.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include "couplingreconstruction.hh"
#include "properties.hh"
#include "plot.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // initialize MPI+X, finalize is done automatically on exit
    Dumux::initialize(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using BulkTypeTag = Properties::TTag::BULKTYPETAG;
    using LowDimTypeTag = Properties::TTag::LOWDIMTYPETAG;

    // try to create a grid (from the given grid file or the input file)
    // for both sub-domains
    using BulkGridManager = Dumux::GridManager<GetPropType<BulkTypeTag, Properties::Grid>>;
    BulkGridManager bulkGridManager;
    bulkGridManager.init("Soil"); // pass parameter group

    using LowDimGridManager = Dumux::GridManager<GetPropType<LowDimTypeTag, Properties::Grid>>;
    LowDimGridManager lowDimGridManager;
    lowDimGridManager.init("Root"); // pass parameter group

    // we compute on the leaf grid view
    const auto& bulkGridView = bulkGridManager.grid().leafGridView();
    const auto& lowDimGridView = lowDimGridManager.grid().leafGridView();

    double hmaxBulk = 0.0, hminBulk = 1e100, hAvgBulk = 0.0;
    std::vector<double> hsBulk(bulkGridView.size(0)); int i = 0;
    for (const auto& element : elements(bulkGridView))
    {
        const auto h = diameter(element.geometry());
        hmaxBulk = std::max(hmaxBulk, h);
        hminBulk = std::min(hminBulk, h);
        hAvgBulk += h;
        hsBulk[i++] = h;
    }

    std::size_t tenPercent = hsBulk.size()/10;
    std::partial_sort(hsBulk.begin(), hsBulk.begin() + tenPercent, hsBulk.end());
    auto hbulkminavg = std::accumulate(hsBulk.begin(), hsBulk.begin() + tenPercent, 0.0);
    hbulkminavg /= double(tenPercent);

    double hmaxLowDim = 0.0, hminLowDim = 1e100, hAvgLowDim = 0.0;
    for (const auto& element : elements(lowDimGridView))
    {
        const auto h = diameter(element.geometry());
        hmaxLowDim = std::max(hmaxLowDim, h);
        hminLowDim = std::min(hminLowDim, h);
        hAvgLowDim += h;
    }

    hAvgBulk /= bulkGridView.size(0);
    hAvgLowDim /= lowDimGridView.size(0);
    std::cout << Fmt::format("Bulk: hmin {:.5e}, hmax {:.5e}, havg {:.5e}\n", hminBulk, hmaxBulk, hAvgBulk);
    std::cout << Fmt::format("Bulk: smallest ten percent havg {:.5e}\n", hbulkminavg);
    std::cout << Fmt::format("Bulk: numElements {}, numVertices {}\n", bulkGridView.size(0), bulkGridView.size(3));
    std::cout << Fmt::format("LowDim: hmin {:.5e}, hmax {:.5e}, havg {:.5e}\n", hminLowDim, hmaxLowDim, hAvgLowDim);
    std::cout << Fmt::format("LowDim: numElements {}, numVertices {}\n", lowDimGridView.size(0), lowDimGridView.size(1));
    if (getParam<bool>("Problem.PrintMeshSizeAndExit", false))
        return 0;

    // create the finite volume grid geometry
    using BulkGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    auto bulkGridGeometry = std::make_shared<BulkGridGeometry>(bulkGridView);
    using LowDimGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
    auto lowDimGridGeometry = std::make_shared<LowDimGridGeometry>(lowDimGridView);

    // the mixed dimension type traits
    using Traits = MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    constexpr auto bulkIdx = Traits::template SubDomain<0>::Index();
    constexpr auto lowDimIdx = Traits::template SubDomain<1>::Index();

    // the coupling manager
    using CouplingManager = GetPropType<BulkTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(bulkGridGeometry, lowDimGridGeometry);

    // the problem (initial and boundary conditions)
    using BulkProblem = GetPropType<BulkTypeTag, Properties::Problem>;
    auto bulkProblem = std::make_shared<BulkProblem>(bulkGridGeometry, couplingManager);

    // the low dim spatial parameters
    using LowDimSpatialParams = GetPropType<LowDimTypeTag, Properties::SpatialParams>;
    auto lowDimSpatialParams = std::make_shared<LowDimSpatialParams>(lowDimGridGeometry, lowDimGridManager.getGridData());

    // the low dim problem (initial and boundary conditions)
    const auto domainSize = bulkGridGeometry->bBoxMax()-bulkGridGeometry->bBoxMin();
    using LowDimProblem = GetPropType<LowDimTypeTag, Properties::Problem>;
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimGridGeometry, lowDimSpatialParams, couplingManager, domainSize);

    // the solution vector
    Traits::SolutionVector sol;
    sol[bulkIdx].resize(bulkGridGeometry->numDofs());
    sol[lowDimIdx].resize(lowDimGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(sol[bulkIdx]);
    lowDimProblem->applyInitialSolution(sol[lowDimIdx]);
    auto oldSol = sol;

    // coupling reconstruction scheme
    using Reconstruction = RootSoil::CouplingReconstruction<typename BulkProblem::SpatialParams::PcKrSwCurve>;
    auto couplingReconstruction = std::make_shared<Reconstruction>();
    couplingManager->init(bulkProblem, lowDimProblem, couplingReconstruction, sol);
    bulkProblem->computePointSourceMap();
    lowDimProblem->computePointSourceMap();

    // enable source cache
    couplingReconstruction->enableCache(*couplingManager);

    // the grid variables
    using BulkGridVariables = GetPropType<BulkTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkGridGeometry);
    bulkGridVariables->init(sol[bulkIdx]);
    using LowDimGridVariables = GetPropType<LowDimTypeTag, Properties::GridVariables>;
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimGridGeometry);
    lowDimGridVariables->init(sol[lowDimIdx]);

    // get some time loop parameters
    using Scalar = Traits::Scalar;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dtInit = getParam<Scalar>("TimeLoop.DtInitial");
    const auto maxTimeStepDivisions = getParam<std::size_t>("Newton.MaxTimeStepDivisions", 10);
    const auto retryTimeStepReductionFactor = getParam<Scalar>("Newton.RetryTimeStepReductionFactor", 0.5);

    const auto enableVtkOutput = getParam<bool>("Output.EnableVtkOutput", true);

    using BulkSolutionVector = std::decay_t<decltype(sol[bulkIdx])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, sol[bulkIdx], bulkProblem->name());
    GetPropType<BulkTypeTag, Properties::IOFields>::initOutputModule(bulkVtkWriter);
    if (enableVtkOutput) bulkVtkWriter.write(0.0);

    using LowDimSolutionVector = std::decay_t<decltype(sol[lowDimIdx])>;
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, sol[lowDimIdx], lowDimProblem->name());
    GetPropType<LowDimTypeTag, Properties::IOFields>::initOutputModule(lowDimVtkWriter);
    std::vector<Scalar> sourceTerms(lowDimGridView.size(0), 0.0);
    std::vector<Scalar> depth(lowDimGridView.size(0), 0.0);
    for (const auto& element : elements(lowDimGridView))
        depth[lowDimGridGeometry->elementMapper().index(element)] = element.geometry().center()[2];

    lowDimVtkWriter.addField(lowDimProblem->spatialParams().getRadii(), "radius");
    lowDimVtkWriter.addField(lowDimProblem->spatialParams().getAges(), "age (days)");
    lowDimVtkWriter.addField(lowDimProblem->spatialParams().getOrders(), "order");
    lowDimVtkWriter.addField(lowDimProblem->spatialParams().getKx(), "Kx (m^4/Pa/s)");
    lowDimVtkWriter.addField(lowDimProblem->spatialParams().getKr(), "Kr (m/Pa/s)");
    lowDimVtkWriter.addField(sourceTerms, "source (kg/s)");
    lowDimVtkWriter.addField(depth, "depth (m)");
    if (enableVtkOutput)
        lowDimVtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dtInit, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(bulkProblem, lowDimProblem),
                                                 std::make_tuple(bulkGridGeometry, lowDimGridGeometry),
                                                 std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                 couplingManager, timeLoop, oldSol);

    // the linear solver
    using LinearSolver = BlockDiagILU0BiCGSTABSolver;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // plotting tools
    RootSoil::TranspirationPlot<LowDimProblem> transpirationPlot(lowDimProblem);
    bool enableTranspirationPlot = getParam<bool>("Problem.PlotTranspiration", true);

    // time loop
    const auto outputInterval = getParam<Scalar>("TimeLoop.EpisodeLength");
    timeLoop->setPeriodicCheckPoint(outputInterval);
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // try solving the non-linear system
        for (std::size_t i = 0; i <= maxTimeStepDivisions; ++i)
        {
            // set time for boundary conditions
            lowDimProblem->setTime(timeLoop->time() + timeLoop->timeStepSize());

            // linearize & solve
            try {
                nonLinearSolver.solve(sol);
                break;
            }
            catch (const NumericalProblem& e)
            {
                if (i < maxTimeStepDivisions)
                {
                    // set solution to previous solution
                    sol = assembler->prevSol();

                    // reset the grid variables to the previous solution
                    assembler->resetTimeStep(sol);

                    std::cout << "Newton solver did not converge with dt = "
                              << timeLoop->timeStepSize() << " seconds. Retrying with time step of "
                              << timeLoop->timeStepSize() * retryTimeStepReductionFactor << " seconds\n";

                    // try again with dt = dt * retryTimeStepReductionFactor
                    timeLoop->setTimeStepSize(timeLoop->timeStepSize() * retryTimeStepReductionFactor);
                }

                else
                {
                    DUNE_THROW(NumericalProblem, "Newton solver didn't converge after "
                                                 << maxTimeStepDivisions << " time-step divisions. dt="
                                                 << timeLoop->timeStepSize() << '\n');
                }
            }
        }

        // make the new solution the old solution
        oldSol = sol;
        bulkGridVariables->advanceTimeStep();
        lowDimGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // output the source terms
        bulkProblem->computeSourceIntegral(sol[bulkIdx], *bulkGridVariables);
        lowDimProblem->computeSourceIntegral(sol[lowDimIdx], *lowDimGridVariables, sourceTerms);

        if (enableTranspirationPlot && (timeLoop->isCheckPoint() || timeLoop->timeStepIndex() == 1))
        {
            const double stepSize = timeLoop->timeStepIndex() == 1 ? dtInit : outputInterval;
            transpirationPlot.addDataPoint(sol[lowDimIdx], *lowDimGridVariables, timeLoop->time(), stepSize);
        }

        // write vtk output
        if (enableVtkOutput && (timeLoop->isCheckPoint() || timeLoop->finished()))
        {
            bulkVtkWriter.write(timeLoop->time());
            lowDimVtkWriter.write(timeLoop->time());
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    }

    timeLoop->finalize(bulkGridGeometry->gridView().comm());

    // print dumux end message
    Parameters::print();

    return 0;

} // end main
