// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Buckley-Leverett test for the two-phase porous-medium flow model.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "analyticsolution.hh"
#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TwoPBuckleyLeverettTpfa;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter);

    BuckleyLeverettAnalyticSolution<TypeTag> analyticSolution(problem);
    vtkWriter.addField(analyticSolution.values(), "Sw_exact");
    vtkWriter.write(0.0);

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    timeLoop->start();
    do {
        nonLinearSolver.solve(x, *timeLoop);

        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        analyticSolution.update(timeLoop->time());
        vtkWriter.write(timeLoop->time());

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    } while (!timeLoop->finished());

    nonLinearSolver.report();
    timeLoop->finalize(leafGridView.comm());

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    constexpr auto saturationIdx = ModelTraits::Indices::saturationIdx;

    Scalar maxSaturationError = 0.0;
    for (std::size_t dofIdx = 0; dofIdx < x.size(); ++dofIdx)
    {
        const auto wettingSaturation = 1.0 - x[dofIdx][saturationIdx];
        maxSaturationError = std::max(maxSaturationError, std::abs(wettingSaturation - analyticSolution.values()[dofIdx]));
    }

    const auto maxAllowedSaturationError = getParam<Scalar>("Problem.MaxSaturationError");
    if (maxSaturationError > maxAllowedSaturationError)
        DUNE_THROW(Dune::InvalidStateException, "Maximum saturation error " << maxSaturationError
                   << " exceeds the threshold " << maxAllowedSaturationError);

    if (leafGridView.comm().rank() == 0)
    {
        std::cout << "Maximum saturation error against Buckley-Leverett solution: "
                  << maxSaturationError << std::endl;
        Parameters::print();
    }

    return 0;
}
