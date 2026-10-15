// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPOneCTests
 * \brief Steady heat pipe with the 2p1cni model compared with the solution of Udell (1985).
 */
#include <config.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include "properties.hh"
#include "heatpipesolution.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TYPETAG;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
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
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(getParam<Scalar>("TimeLoop.MaxTimeStepSize"));

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // the steady state is reached by implicit Euler steps of growing size
    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);
        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();
        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    } while (!timeLoop->finished());
    timeLoop->finalize(leafGridView.comm());

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(timeLoop->time());

    // discrete L2 errors with the control volumes as weights
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    const TwoPOneCHeatPipeSolution<Scalar> solution;
    Scalar errorS = 0.0, errorPg = 0.0, errorPl = 0.0, errorT = 0.0;
    Scalar meshSize = 0.0;
    std::ofstream profile(problem->name() + "-profile.csv");
    profile << "x,saturation,gas_pressure,liquid_pressure,temperature,"
            << "saturation_exact,gas_pressure_exact,liquid_pressure_exact,temperature_exact\n"
            << std::setprecision(15);

    auto fvGeometry = localView(*gridGeometry);
    auto elemVolVars = localView(gridVariables->curGridVolVars());
    for (const auto& element : elements(leafGridView))
    {
        using std::max;
        meshSize = max(meshSize, element.geometry().volume());
        fvGeometry.bindElement(element);
        elemVolVars.bindElement(element, fvGeometry, x);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            const Scalar pos = scv.dofPosition()[0];
            const auto exact = solution(pos);
            const Scalar s = volVars.saturation(FluidSystem::liquidPhaseIdx);
            const Scalar pg = volVars.pressure(FluidSystem::gasPhaseIdx);
            const Scalar pl = volVars.pressure(FluidSystem::liquidPhaseIdx);
            const Scalar T = volVars.temperature();

            errorS += (s - exact.saturation)*(s - exact.saturation)*scv.volume();
            errorPg += (pg - exact.gasPressure)*(pg - exact.gasPressure)*scv.volume();
            errorPl += (pl - exact.liquidPressure)*(pl - exact.liquidPressure)*scv.volume();
            errorT += (T - exact.temperature)*(T - exact.temperature)*scv.volume();

            profile << pos << "," << s << "," << pg << "," << pl << "," << T << ","
                    << exact.saturation << "," << exact.gasPressure << ","
                    << exact.liquidPressure << "," << exact.temperature << "\n";
        }
    }

    using std::sqrt;
    std::cout << std::setprecision(6) << std::scientific
              << "L2 errors for mesh size " << meshSize << " m: saturation " << sqrt(errorS)
              << ", gas pressure " << sqrt(errorPg) << " Pa, liquid pressure " << sqrt(errorPl)
              << " Pa, temperature " << sqrt(errorT) << " K" << std::endl;

    std::ofstream errors(problem->name() + "-errors.csv", std::ios::app);
    errors << std::setprecision(15) << meshSize << ","
           << sqrt(errorS) << "," << sqrt(errorPg) << "," << sqrt(errorPl) << "," << sqrt(errorT) << "\n";

    return 0;
}
