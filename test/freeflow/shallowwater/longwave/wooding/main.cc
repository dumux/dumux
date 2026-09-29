// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Rainfall runoff on Wooding's V-catchment
 *
 * Checks that the water balance of the catchment closes, and optionally that the outflow at
 * the end of the rain equals the rainfall on the catchment, i.e. that the catchment has
 * reached equilibrium.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include <dune/common/exceptions.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/timeserieswriter.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::TYPETAG;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    GridManager<Grid> gridManager;
    gridManager.init();
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());

    auto problem = std::make_shared<Problem>(gridGeometry);

    SolutionVector sol;
    problem->applyInitialSolution(sol);
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    IOFields::initOutputModule(vtkWriter);

    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(
        0.0, getParam<Scalar>("TimeLoop.DtInitial"), getParam<Scalar>("TimeLoop.TEnd")
    );
    timeLoop->setMaxTimeStepSize(getParam<Scalar>("TimeLoop.MaxTimeStepSize"));

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = UMFPackIstlSolver<LinearSolverTraits<GridGeometry>,
                                           LinearAlgebraTraitsFromAssembler<Assembler>>;
    NewtonSolver<Assembler, LinearSolver> solver(
        std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, sol),
        std::make_shared<LinearSolver>()
    );

    const auto forEachScv = [&](const auto& curSol, auto&& f)
    {
        auto fvGeometry = localView(*gridGeometry);
        auto elemVolVars = localView(gridVariables->curGridVolVars());
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, curSol);
            for (const auto& scv : scvs(fvGeometry))
                f(scv, elemVolVars[scv]);
        }
    };

    const auto storage = [&](const auto& curSol)
    {
        Scalar volume = 0.0;
        forEachScv(curSol, [&](const auto& scv, const auto& volVars){ volume += volVars.waterDepth()*scv.volume(); });
        return volume;
    };

    Scalar catchmentArea = 0.0;
    forEachScv(sol, [&](const auto& scv, const auto&){ catchmentArea += scv.volume(); });

    const auto outletDischarge = [&](const auto& curSol)
    {
        Scalar discharge = 0.0;
        auto fvGeometry = localView(*gridGeometry);
        auto elemVolVars = localView(gridVariables->curGridVolVars());
        auto elemFluxVarsCache = localView(gridVariables->gridFluxVarsCache());
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            if (!element.hasBoundaryIntersections())
                continue;

            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, curSol);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);
            for (const auto& scvf : scvfs(fvGeometry))
                if (scvf.boundary())
                    discharge += problem->neumann(element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf)[Indices::massBalanceIdx]
                                 *scvf.area();
        }
        return discharge;
    };

    // mean channel depth over a band of degrees of freedom near the outlet
    const auto probeMinY = getParam<Scalar>("Problem.DepthProbeMinY");
    const auto probeMaxY = getParam<Scalar>("Problem.DepthProbeMaxY");
    const auto outletDepth = [&](const auto& curSol)
    {
        Scalar depth = 0.0, volume = 0.0;
        forEachScv(curSol, [&](const auto& scv, const auto& volVars)
        {
            const auto& pos = scv.dofPosition();
            if (pos[1] < probeMinY || pos[1] > probeMaxY || !problem->inChannel(pos))
                return;
            depth += volVars.waterDepth()*scv.volume();
            volume += scv.volume();
        });
        return volume > 0.0 ? depth/volume : 0.0;
    };

    TimeSeriesWriter hydrograph(problem->name() + "_hydrograph.dat", "discharge", "m^3/s");
    TimeSeriesWriter stage(problem->name() + "_depth.dat", "depth", "m");
    hydrograph.write(0.0, outletDischarge(sol));
    stage.write(0.0, outletDepth(sol));

    const auto rainEventEnd = getParam<Scalar>("Problem.RainEventEnd");
    const auto rainFallRate = getParam<Scalar>("Problem.RainFallRate");
    timeLoop->setCheckPoint(rainEventEnd);
    const auto initialStorage = storage(sol);
    Scalar totalRain = 0.0, totalOutflow = 0.0, dischargeAtRainEnd = 0.0;

    vtkWriter.write(0.0);
    timeLoop->start(); do
    {
        const bool raining = timeLoop->time() < rainEventEnd;
        problem->setRainActive(raining);
        solver.solve(sol, *timeLoop);

        const auto dt = timeLoop->timeStepSize();
        const auto discharge = outletDischarge(sol);
        totalRain += problem->rainFallRate()*catchmentArea*dt;
        totalOutflow += discharge*dt;
        if (raining)
            dischargeAtRainEnd = discharge;

        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        hydrograph.write(timeLoop->time(), discharge);
        stage.write(timeLoop->time(), outletDepth(sol));

        if (timeLoop->finished())
            vtkWriter.write(timeLoop->time());

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(solver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(gridGeometry->gridView().comm());

    using std::abs;
    const auto balanceError = abs(totalRain - (storage(sol) - initialStorage) - totalOutflow)/totalRain;
    std::cout << "Water balance: rain " << totalRain << " m^3, outflow " << totalOutflow
              << " m^3, relative error " << balanceError << std::endl;

    const auto balanceTolerance = getParam<Scalar>("Problem.WaterBalanceTolerance");
    if (!(balanceError <= balanceTolerance))
        DUNE_THROW(Dune::Exception, "Water balance not closed: relative error " << balanceError
                                    << " exceeds " << balanceTolerance);

    if (hasParam("Problem.EquilibriumTolerance"))
    {
        const auto equilibriumDischarge = rainFallRate*catchmentArea;
        const auto equilibriumError = abs(dischargeAtRainEnd/equilibriumDischarge - 1.0);
        std::cout << "Discharge at the end of the rain: " << dischargeAtRainEnd << " m^3/s, equilibrium "
                  << equilibriumDischarge << " m^3/s, relative deviation " << equilibriumError << std::endl;

        const auto equilibriumTolerance = getParam<Scalar>("Problem.EquilibriumTolerance");
        if (!(equilibriumError <= equilibriumTolerance))
            DUNE_THROW(Dune::Exception, "Equilibrium not reached: relative deviation " << equilibriumError
                                        << " exceeds " << equilibriumTolerance);
    }

    return 0;
}
