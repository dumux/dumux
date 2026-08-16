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
#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

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
#include <dumux/common/integrate.hh>

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
    SolutionVector sol(gridGeometry->numDofs());
    problem->applyInitialSolution(sol);
    auto solOld = sol;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
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
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, solOld);

    using LinearSolver = AMGBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    timeLoop->start();
    do {
        nonLinearSolver.solve(sol, *timeLoop);

        solOld = sol;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        analyticSolution.update(timeLoop->time());
        vtkWriter.write(timeLoop->time());

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    } while (!timeLoop->finished());

    nonLinearSolver.report();
    timeLoop->finalize(leafGridView.comm());


    // compute relative error in wetting-phase center of mass and total mass
    // assumptions: solutions are (pseudo)1D, densities and porosity are constant, TPFA discretization
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    constexpr auto saturationIdx = ModelTraits::Indices::saturationIdx;

    FluidState fluidState;
    Scalar referencePressure = getParam<Scalar>("Problem.ReferencePressure");
    fluidState.setTemperature(problem->spatialParams().temperatureAtPos(GlobalPosition{}));
    fluidState.setPressure(FluidSystem::phase0Idx, referencePressure);
    fluidState.setPressure(FluidSystem::phase1Idx, referencePressure);
    const auto densityW = FluidSystem::density(fluidState, FluidSystem::phase0Idx);
    const auto porosity = problem->spatialParams().porosityAtPos(GlobalPosition{});
    const auto xMin = gridGeometry->bBoxMin()[0];
    const auto xMax = gridGeometry->bBoxMax()[0];
    const auto domainWidth = gridGeometry->bBoxMax()[1] - gridGeometry->bBoxMin()[1];

    // analytic computation of center of mass and total mass
    auto swFunc = [&analyticSolution, &tEnd](Scalar x)
    {
        return analyticSolution.computeSaturation(x, tEnd);
    };

    auto firstMomentOfMassWIntegrand = [&densityW, &porosity, &domainWidth, &swFunc](Scalar x)
    {
        return x*densityW*swFunc(x)*porosity*domainWidth;
    };
    Scalar firstMomentOfMassWAnalytic = integrateScalarFunction(firstMomentOfMassWIntegrand, xMin, xMax);

    auto totalMassWIntegrand = [&densityW, &porosity, &domainWidth, &swFunc](Scalar x)
    {
        return densityW*swFunc(x)*porosity*domainWidth;
    };
    Scalar totalMassWAnalytic = integrateScalarFunction(totalMassWIntegrand, xMin, xMax);

    // numeric computation of center of mass and total mass
    Scalar firstMomentOfMassWNumeric = 0.0;
    Scalar totalMassWNumeric = 0.0;
    for (const auto& element : elements(gridGeometry->gridView(), Dune::Partitions::interior))
    {
        const auto globalPos = element.geometry().center();
        const auto eIdx = gridGeometry->elementMapper().index(element);
        const auto volume = element.geometry().volume();

        Scalar satWNumeric = 1.0-sol[eIdx][saturationIdx];
        Scalar localMassWNumeric = densityW * porosity * volume * satWNumeric;

        firstMomentOfMassWNumeric += globalPos[0] * localMassWNumeric;
        totalMassWNumeric += localMassWNumeric;
    }

    if (leafGridView.comm().size() > 1)
    {
        firstMomentOfMassWNumeric = leafGridView.comm().sum(firstMomentOfMassWNumeric);
        totalMassWNumeric = leafGridView.comm().sum(totalMassWNumeric);
    }

    Scalar centerOfMassNumeric = firstMomentOfMassWNumeric/totalMassWNumeric;
    Scalar centerOfMassAnalytic = firstMomentOfMassWAnalytic/totalMassWAnalytic;

    // compute relative errors
    Scalar distanceCenterOfMass = std::abs(centerOfMassAnalytic - centerOfMassNumeric);
    Scalar relErrorCenterOfMass = distanceCenterOfMass/centerOfMassAnalytic;
    Scalar differenceTotalMass = std::abs(totalMassWAnalytic - totalMassWNumeric);
    Scalar relErrorTotalMass = differenceTotalMass/totalMassWAnalytic;
    const auto maxError = getParam<Scalar>("Problem.MaxRelError");
    const bool centerOfMassFailed = relErrorCenterOfMass > maxError;
    const bool totalMassFailed = relErrorTotalMass > maxError;
    if (centerOfMassFailed || totalMassFailed)
    {
        std::ostringstream msg;
        msg << "Analytical Buckley-Leverett check failed.";

        if (centerOfMassFailed)
            msg << " Relative wetting-phase center of mass error "
                << relErrorCenterOfMass << " exceeds the threshold " << maxError << ".";

        if (totalMassFailed)
            msg << " Relative total wetting-phase mass error "
                << relErrorTotalMass << " exceeds the threshold " << maxError << ".";

        DUNE_THROW(Dune::InvalidStateException, msg.str());
    }

    if (leafGridView.comm().rank() == 0)
    {
        std::cout << "numeric center of mass: " << centerOfMassNumeric << std::endl;
        std::cout << "analytic center of mass: " << centerOfMassAnalytic << std::endl;
        std::cout << "integrated numeric mass: " << totalMassWNumeric << std::endl;
        std::cout << "integrated analytics mass: " << totalMassWAnalytic << std::endl;
        Parameters::print();
    }

    return 0;
}
