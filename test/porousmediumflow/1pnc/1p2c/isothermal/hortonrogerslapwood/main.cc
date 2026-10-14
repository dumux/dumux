// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePNCTests
 * \brief Growth rate of a convection cell in the Horton-Rogers-Lapwood problem compared with linear stability theory
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::HortonRogersLapwood;

    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    constexpr int soluteIdx = FluidSystem::comp1Idx;

    // Rayleigh number Ra = K drho g H / (mu D_eff) and time scale phi H^2 / D_eff of the homogeneous layer
    const Scalar height = gridGeometry->bBoxMax()[1] - gridGeometry->bBoxMin()[1];
    const Scalar length = gridGeometry->bBoxMax()[0] - gridGeometry->bBoxMin()[0];
    auto fvGeometry = localView(*gridGeometry);
    auto elemVolVars = localView(gridVariables->curGridVolVars());
    const auto element = *elements(leafGridView).begin();
    fvGeometry.bind(element);
    elemVolVars.bind(element, fvGeometry, x);
    const auto& volVars = elemVolVars[*scvs(fvGeometry).begin()];
    const Scalar effectiveDiffusivity = volVars.effectiveDiffusionCoefficient(0, FluidSystem::comp0Idx, soluteIdx);
    const Scalar g = -problem->spatialParams().gravity(gridGeometry->bBoxMin())[1];
    const Scalar rayleigh = volVars.permeability()*problem->densityDifference()*g*height
                            /(volVars.viscosity(0)*effectiveDiffusivity);
    const Scalar timeScale = volVars.porosity()*height*height/effectiveDiffusivity;

    // linear stability theory: sigma = Ra k^2/(k^2 + pi^2) - (k^2 + pi^2), k in units of 1/H, sigma in units of 1/timeScale
    const Scalar k = M_PI*height/length;
    const Scalar diffusiveRate = k*k + M_PI*M_PI;
    const Scalar exactGrowthRate = rayleigh*k*k/diffusiveRate - diffusiveRate;

    // projection of the deviation from the linear profile onto the normal mode
    const auto modeAmplitudeAndMaxDeviation = [&]()
    {
        Scalar num = 0.0, den = 0.0, maxDeviation = 0.0;
        for (const auto& element : elements(leafGridView))
        {
            fvGeometry.bindElement(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& pos = scv.dofPosition();
                const Scalar deviation = x[scv.dofIndex()][soluteIdx] - problem->baseProfile(pos);
                const Scalar mode = problem->mode(pos);
                num += deviation*mode*scv.volume();
                den += mode*mode*scv.volume();
                maxDeviation = std::max(maxDeviation, std::abs(deviation));
            }
        }
        return std::make_pair(num/den, maxDeviation);
    };

    const Scalar dt = getParam<Scalar>("TimeLoop.DimensionlessTimeStepSize")*timeScale;
    const auto numTimeSteps = getParam<int>("TimeLoop.NumTimeSteps");
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, numTimeSteps*dt);
    timeLoop->setMaxTimeStepSize(dt);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    VtkOutputModule vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    std::ofstream logFile(problem->name() + "_growthrate.csv");
    logFile << std::setprecision(12) << "time,amplitude,growthRate,maxDeviation\n";
    auto [amplitudeOld, maxDeviation] = modeAmplitudeAndMaxDeviation();
    logFile << 0.0 << "," << amplitudeOld << ",," << maxDeviation << "\n";

    Scalar growthRate = 0.0;
    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);
        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        const auto [amplitude, deviation] = modeAmplitudeAndMaxDeviation();
        maxDeviation = deviation;
        // an implicit Euler step multiplies a mode with growth rate sigma by 1/(1 - sigma dt)
        if (amplitude != 0.0)
            growthRate = (1.0 - amplitudeOld/amplitude)/dt*timeScale;
        logFile << timeLoop->time()/timeScale << "," << amplitude << "," << growthRate << "," << maxDeviation << "\n";
        amplitudeOld = amplitude;

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(dt);
    } while (!timeLoop->finished());

    vtkWriter.write(timeLoop->time());
    timeLoop->finalize(leafGridView.comm());

    const Scalar error = std::abs(growthRate - exactGrowthRate)/diffusiveRate;
    std::cout << std::setprecision(10)
              << "Ra = " << rayleigh << ", k = " << k
              << ", growth rate " << growthRate << ", exact " << exactGrowthRate
              << ", error " << error << ", max deviation from linear profile " << maxDeviation << std::endl;

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    if (hasParam("Problem.GrowthRateTolerance") && !(error < getParam<Scalar>("Problem.GrowthRateTolerance")))
        DUNE_THROW(Dune::Exception, "Growth rate error " << error << " exceeds the tolerance");

    return 0;
}
