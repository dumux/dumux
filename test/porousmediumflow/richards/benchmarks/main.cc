// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RichardsTests
 * \brief Richards benchmark test
 *
 * Infiltration benchmark:
 * Root-soil benchmark paper Schnepf et al. (case M2.1, Eq. 4) https://doi.org/10.3389/fpls.2020.00316
 * based on Vanderborght 2005 (see Fig. 4abc and Eq. 56-60) https://doi.org/10.2113/4.1.206
 *
 * Evaporation benchmark:
 * Root-soil benchmark paper Schnepf et al. (case M2.2) https://doi.org/10.3389/fpls.2020.00316
 * based on Vanderborght 2005 (see Fig. 5abcd and Eq. 39-47) https://doi.org/10.2113/4.1.206
 */

#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/math.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>

#include <dumux/io/format.hh>
#include <dumux/io/gnuplotinterface.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/porousmediumflow/richards/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <test/porousmediumflow/richards/benchmarks/analytical.hh>
#include <test/porousmediumflow/richards/benchmarks/properties.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Choose tpfa scheme
    using TypeTag = Properties::TTag::RichardsBenchmarkCCTpfa;

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // initialize the vtk output module
    auto vtkWriter = std::make_shared<VtkOutputModule<GridVariables, SolutionVector>>(*gridVariables, x, problem->name());
    vtkWriter->addVolumeVariable([](const auto& volVars){ return volVars.saturation(0)*volVars.porosity(); }, "water content");
    vtkWriter->addVolumeVariable([](const auto& volVars){ return volVars.saturation(0); }, "saturation");
    vtkWriter->addVolumeVariable([](const auto& volVars){ return volVars.pressure(0); }, "pressure");

    const auto dt = getParam<double>("TimeLoop.DtInitial");
    const auto checkPoints = getParam<std::vector<double>>("TimeLoop.TEnd");
    const auto tEnd = checkPoints.back();

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<double>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(getParam<double>("TimeLoop.MaxTimeStepSize"));
    timeLoop->setCheckPoint(checkPoints.begin(), checkPoints.end());
    int checkPointCounter = 0;

    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    using Newton = RichardsNewtonSolver<Assembler, LinearSolver>;
    Newton nonLinearSolver(assembler, linearSolver);

    // find out which scenario we are running (evaporation/infiltration)
    const auto soilType = getParam<std::string>("SpatialParams.Name");
    const double rate = getParam<double>("Problem.SurfaceFluxMilliMeterPerDay");
    const BenchmarkScenario scenario = (rate > 0) ? BenchmarkScenario::evaporation : BenchmarkScenario::infiltration;

    // compute analytical solution for given scenario
    std::unique_ptr<AnalyticalSolutionM21> solutionInfiltration;
    std::unique_ptr<AnalyticalSolutionM22> solutionEvaporation;
    if (scenario == BenchmarkScenario::infiltration)
        solutionInfiltration = std::make_unique<AnalyticalSolutionM21>();
    else if (scenario == BenchmarkScenario::evaporation)
        solutionEvaporation = std::make_unique<AnalyticalSolutionM22>();

    ///////////////////////////////////////
    // for evaporation scenario output
    ///////////////////////////////////////

    const bool openGnuplotWindow = getParam<bool>("Problem.OpenGnuplotWindow", false);
    GnuplotInterface<double> gnuplotEvap(openGnuplotWindow);
    gnuplotEvap.setOpenPlotWindow(openGnuplotWindow);
    std::vector<double> outTime, outActual, outAnalytic;

    ///////////////////////////////////////
    ///////////////////////////////////////

    ///////////////////////////////////////
    // for infiltration scenario output
    ///////////////////////////////////////

    const auto updateAnalyticalWaterContent = [&](auto& theta)
    {
        const auto t = timeLoop->time()/86400.0;
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = leafGridView.indexSet().index(element);
            const auto pos = element.geometry().center()[GridGeometry::GridView::dimensionworld -1];
            theta[eIdx] = solutionInfiltration->waterContent(pos*100, t);
        }
    };

    const auto updateDeltaEtaOutput = [&](auto& deltaEtaNum, auto& thetaNum)
    {
        const auto t = timeLoop->time()/86400.0;
        auto fvGeometry = localView(*gridGeometry);
        auto elemVolVars = localView(gridVariables->curGridVolVars());
        for (const auto& element : elements(leafGridView))
        {
            const auto eIdx = leafGridView.indexSet().index(element);
            const auto pos = element.geometry().center()[GridGeometry::GridView::dimensionworld -1];
            deltaEtaNum[eIdx] = solutionInfiltration->deltaEta(pos*100, t);

            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, x);
            for (const auto& scv : scvs(fvGeometry))
                thetaNum[eIdx] = elemVolVars[scv].saturation() * elemVolVars[scv].porosity();
        }
    };

    std::vector<double> pos(leafGridView.size(0));
    for (const auto& element : elements(leafGridView))
        pos[leafGridView.indexSet().index(element)] = element.geometry().center()[GridGeometry::GridView::dimensionworld -1];


    GnuplotInterface<double> gnuplotInfilt(openGnuplotWindow);
    gnuplotInfilt.setOpenPlotWindow(openGnuplotWindow);
    std::vector<double> deltaEtaNum(leafGridView.size(0)), thetaNum(leafGridView.size(0));

    std::vector<double> analyticalWaterContent;
    if (scenario == BenchmarkScenario::infiltration)
    {
        analyticalWaterContent.resize(leafGridView.size(0), solutionInfiltration->initialWaterContent());
        vtkWriter->addField(analyticalWaterContent, "analytic water content");

        // analytical solution
        const auto& [theta, deltaEta] = solutionInfiltration->waterContentAndDeltaEtaPlot();
        gnuplotInfilt.addDataSetToPlot(theta, deltaEta,
            Fmt::format("theta_deltaeta_ana_{}_{:d}.dat", soilType, checkPointCounter), "w l t 'analytic'"
        );

        const auto yLimits = getParam<std::vector<double>>("Analytical.PlotLimits");
        gnuplotInfilt.setYRange(yLimits[0], yLimits[1]);
        gnuplotInfilt.setXRange(0.0, 0.5);
    }

    ///////////////////////////////////////
    ///////////////////////////////////////

    // initial solution
    vtkWriter->write(0.0);

    // time loop
    timeLoop->start();
    while (!timeLoop->finished())
    {
        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        ///////////////////////////////////////
        // output for infiltration scenario
        ///////////////////////////////////////
        if (scenario == BenchmarkScenario::infiltration && (timeLoop->finished() || timeLoop->isCheckPoint()))
        {
            ++checkPointCounter;

            // write analytical solution to ParaView
            updateAnalyticalWaterContent(analyticalWaterContent);
            vtkWriter->write(timeLoop->time());

            updateDeltaEtaOutput(deltaEtaNum, thetaNum);

            const auto posRef = interpolate<InterpolationPolicy::LinearTable>(solutionInfiltration->refWaterContent(), thetaNum, pos);
            std::cout << Fmt::format("Simulated x_a: {}cm\n", posRef);

            gnuplotInfilt.addDataSetToPlot(thetaNum, deltaEtaNum,
                Fmt::format("theta_deltaeta_num_{}_{:d}.dat", soilType, checkPointCounter),
                Fmt::format("w p t 'DuMux {:.1f} days'", timeLoop->time()/86400.0)
            );

            gnuplotInfilt.plot(Fmt::format("infiltration_{}_{:d}", soilType, checkPointCounter));
        }

        ///////////////////////////////////////
        // output for evaporation scenario
        ///////////////////////////////////////
        if (scenario == BenchmarkScenario::evaporation)
        {
            if (timeLoop->finished())
                vtkWriter->write(timeLoop->time());

            const double rate = problem->computeActualRate(x, *gridVariables, false);
            const double analyticRate = solutionEvaporation->evaporationRate(timeLoop->time()/86400.0);
            std::cout << Fmt::format("Actual rate: {:.5g} (mm/day), analytic rate: {:.5g} (mm/day)\n", rate, analyticRate);

            outTime.push_back(timeLoop->time()/86400.0);
            outActual.push_back(rate);
            outAnalytic.push_back(analyticRate);

            gnuplotEvap.resetPlot();
            gnuplotEvap.setXRange(0, tEnd/86400.0);
            gnuplotEvap.setXlabel("time in days");
            gnuplotEvap.setYlabel("evaporation rate in mm/days");
            gnuplotEvap.addDataSetToPlot(outTime, outActual, Fmt::format("rate_actual_{}.dat", soilType), "w p t 'DuMux'");
            gnuplotEvap.addDataSetToPlot(outTime, outAnalytic, Fmt::format("rate_analytical_{}.dat", soilType), "w l t 'reference'");
            gnuplotEvap.plot(Fmt::format("evaporation_{}", soilType));
        }

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    }

    timeLoop->finalize(leafGridView.comm());

    nonLinearSolver.report();

    // print dumux end message
    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;
}
