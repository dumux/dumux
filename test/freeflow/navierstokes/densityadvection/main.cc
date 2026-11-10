// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the advective momentum flux of the CVFE momentum model with variable density.
 *
 * A density interface is translated with a uniform velocity. The momentum balance with a constant
 * pressure and without viscosity has the uniform velocity as solution. The maximum deviation of the
 * discrete velocity from it is reported for several density ratios and grid resolutions.
 */

#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/assembler.hh>

#include "properties.hh"

#ifndef TYPETAG_DENSITY_ADVECTION
#define TYPETAG_DENSITY_ADVECTION DensityAdvectionTestPQ1BubbleHybrid
#endif

namespace Dumux {

using TypeTag = Properties::TTag::TYPETAG_DENSITY_ADVECTION;
using Scalar = GetPropType<TypeTag, Properties::Scalar>;
using GridType = GetPropType<TypeTag, Properties::Grid>;
using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
using Problem = GetPropType<TypeTag, Properties::Problem>;
using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

//! Maximum deviation of the velocity from the uniform velocity at the end time
Scalar maxVelocityDeviation(int numCellsX, Scalar baseDensity, Scalar densityRatio, Scalar cfl, const std::string& densityField)
{
    const auto length = getParam<Scalar>("Grid.Length");
    const auto height = getParam<Scalar>("Grid.Height");
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const int numCellsY = std::max(1, static_cast<int>(std::round(numCellsX*height/length)));

    GridType grid({0.0, 0.0}, {length, height}, {numCellsX, numCellsY});
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());

    auto problem = std::make_shared<Problem>(gridGeometry);
    problem->setDensities(baseDensity, baseDensity*densityRatio);
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    const Scalar velocity = problem->exactVelocity()[0];
    const Scalar h = length/numCellsX;
    const int numSteps = static_cast<int>(std::ceil(tEnd*velocity/(cfl*h)));
    const Scalar dt = tEnd/numSteps;
    if (densityField != "exact")
        problem->computeTransportedDensity(densityField == "upwind", numCellsX, length, dt, numSteps);
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd, false);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    timeLoop->start(); do
    {
        problem->setTimeStep(timeLoop->time(), timeLoop->time() + timeLoop->timeStepSize());
        nonLinearSolver.solve(x);
        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();
    } while (!timeLoop->finished());

    // only the control-volume dofs carry nodal velocities, hybrid dofs hold corrections
    Scalar maxDeviation = 0.0;
    auto elemDisc = localView(*gridGeometry);
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        elemDisc.bind(element);
        for (const auto& scv : scvs(elemDisc))
        {
            using std::abs;
            using std::isfinite;
            const auto& values = x[scv.dofIndex()];
            if (!isfinite(values[0]) || !isfinite(values[1]))
                return std::numeric_limits<Scalar>::quiet_NaN();
            maxDeviation = std::max({maxDeviation, abs(values[0] - velocity), abs(values[1])});
        }
    }
    return maxDeviation;
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    const auto baseCells = getParam<int>("Study.BaseCells");
    const auto numLevels = getParam<int>("Study.NumLevels");
    const auto cfl = getParam<Scalar>("Study.Cfl");
    const auto baseDensity = getParam<Scalar>("Study.BaseDensity", 1.0);
    const auto ratios = getParam<std::vector<Scalar>>("Study.DensityRatios");
    const auto densityField = getParam<std::string>("Study.DensityField", "exact");

    // optional assertions on the finest grid, one entry per density ratio
    const auto maxErrors = getParam<std::vector<Scalar>>("Study.MaxError", std::vector<Scalar>{});
    const auto minOrders = getParam<std::vector<Scalar>>("Study.MinOrder", std::vector<Scalar>{});

    bool passed = true;
    std::cout << std::setprecision(6) << std::scientific;
    for (std::size_t ratioIdx = 0; ratioIdx < ratios.size(); ++ratioIdx)
    {
        const auto ratio = ratios[ratioIdx];
        Scalar previous = 0.0;
        Scalar order = 0.0;
        Scalar deviation = 0.0;
        for (int level = 0; level < numLevels; ++level)
        {
            const int numCellsX = baseCells << level;
            deviation = maxVelocityDeviation(numCellsX, baseDensity, ratio, cfl, densityField);
            std::cout << "[Result] density " << densityField << " ratio " << ratio << " cells " << numCellsX
                      << " max|u-u0| " << deviation;
            if (level > 0 && deviation > 0 && previous > 0)
            {
                order = std::log2(previous/deviation);
                std::cout << " order " << order;
            }
            std::cout << std::endl;
            previous = deviation;
        }

        if (ratioIdx < maxErrors.size() && !(deviation <= maxErrors[ratioIdx]))
        {
            std::cout << "[Failed] max|u-u0| " << deviation << " exceeds " << maxErrors[ratioIdx] << std::endl;
            passed = false;
        }
        if (ratioIdx < minOrders.size() && !(order >= minOrders[ratioIdx]))
        {
            std::cout << "[Failed] observed order " << order << " is below " << minOrders[ratioIdx] << std::endl;
            passed = false;
        }
    }

    return passed ? 0 : 1;
}
