// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief McWhorter-Sunada test for the two-phase porous-medium flow model.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <sstream>

#include <dune/common/exceptions.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/integrate.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"
#include "analyticsolution.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TwoPMcWhorterSunadaTpfa;

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
    if (!(tEnd > 0.0) || !(dt > 0.0) || !(maxDt > 0.0)
        || !(getParam<Scalar>("Problem.MaxRelError") > 0.0))
        DUNE_THROW(Dune::InvalidStateException, "McWhorter benchmark requires positive times and error tolerance");

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector sol(gridGeometry->numDofs());
    problem->applyInitialSolution(sol);
    auto solOld = sol;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    McWhorterAnalyticSolution<TypeTag> reference(problem);
    if (reference.computeSaturation(gridGeometry->bBoxMax()[0], tEnd) > reference.initialSaturation())
        DUNE_THROW(Dune::InvalidStateException,
                   "McWhorter reference front reaches the closed right boundary at TEnd. "
                   "Reduce TEnd or enlarge the domain.");
    reference.update(0.0);
    std::vector<Scalar> saturationError(gridGeometry->numDofs(), 0.0);

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.addField(reference.values(), "S_aq_reference");
    vtkWriter.addField(saturationError, "S_aq_error");
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

        reference.update(timeLoop->time());
        constexpr auto saturationIdx = GetPropType<TypeTag, Properties::ModelTraits>::Indices::saturationIdx;
        for (std::size_t i = 0; i < sol.size(); ++i)
            saturationError[i] = (1.0 - sol[i][saturationIdx]) - reference.values()[i];
        vtkWriter.write(timeLoop->time());

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    } while (!timeLoop->finished());

    nonLinearSolver.report();
    timeLoop->finalize(leafGridView.comm());

    // Compute relative errors in imbibed wetting-phase mass, its center of mass,
    // and the saturation profile. The comparison assumes a pseudo-1D domain with
    // constant density and porosity, as required by the semi-analytical reference.
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    constexpr auto saturationIdx = ModelTraits::Indices::saturationIdx;

    FluidState fluidState;
    const Scalar referencePressure = problem->referencePressure();
    fluidState.setTemperature(problem->spatialParams().temperatureAtPos(GlobalPosition{}));
    fluidState.setPressure(FluidSystem::phase0Idx, referencePressure);
    fluidState.setPressure(FluidSystem::phase1Idx, referencePressure);
    const Scalar densityW = FluidSystem::density(fluidState, FluidSystem::phase0Idx);
    const Scalar porosity = problem->spatialParams().porosityAtPos(GlobalPosition{});
    const Scalar swr = reference.initialSaturation();
    const Scalar xMin = gridGeometry->bBoxMin()[0];
    const Scalar xMax = gridGeometry->bBoxMax()[0];
    const Scalar domainWidth = gridGeometry->bBoxMax()[1] - gridGeometry->bBoxMin()[1];

    // Integrate the semi-analytical reference over the finite computational domain.
    auto swReference = [&reference, &swr, &tEnd](Scalar x)
    {
        return reference.computeSaturation(x, tEnd) - swr;
    };
    auto imbibedMassIntegrand = [&densityW, &porosity, &domainWidth, &swReference](Scalar x)
    {
        return densityW*porosity*domainWidth*swReference(x);
    };
    auto firstMomentIntegrand = [&xMin, &imbibedMassIntegrand](Scalar x)
    {
        return (x - xMin)*imbibedMassIntegrand(x);
    };
    const Scalar imbibedMassReference = integrateScalarFunction(imbibedMassIntegrand, xMin, xMax);
    const Scalar firstMomentReference = integrateScalarFunction(firstMomentIntegrand, xMin, xMax);

    // Integrate the numerical solution over interior cells.
    Scalar imbibedMassNumeric = 0.0;
    Scalar firstMomentNumeric = 0.0;
    Scalar saturationL1Error = 0.0;
    for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
    {
        const auto index = gridGeometry->elementMapper().index(element);
        const Scalar swNumeric = 1.0 - sol[index][saturationIdx];
        const Scalar imbibedMass = densityW*porosity*element.geometry().volume()*(swNumeric - swr);
        imbibedMassNumeric += imbibedMass;
        firstMomentNumeric += (element.geometry().center()[0] - xMin)*imbibedMass;
        saturationL1Error += densityW*porosity*element.geometry().volume()
                             *std::abs(swNumeric - reference.values()[index]);
    }

    if (leafGridView.comm().size() > 1)
    {
        imbibedMassNumeric = leafGridView.comm().sum(imbibedMassNumeric);
        firstMomentNumeric = leafGridView.comm().sum(firstMomentNumeric);
        saturationL1Error = leafGridView.comm().sum(saturationL1Error);
    }

    const Scalar centerOfMassNumeric = firstMomentNumeric/imbibedMassNumeric;
    const Scalar centerOfMassReference = firstMomentReference/imbibedMassReference;
    const Scalar relativeCenterOfMassError = std::abs(centerOfMassReference - centerOfMassNumeric)
                                             /centerOfMassReference;
    const Scalar relativeImbibedMassError = std::abs(imbibedMassReference - imbibedMassNumeric)
                                            /imbibedMassReference;
    const Scalar relativeSaturationL1Error = saturationL1Error/imbibedMassReference;
    const Scalar maxError = getParam<Scalar>("Problem.MaxRelError");
    const bool centerOfMassFailed = relativeCenterOfMassError > maxError;
    const bool imbibedMassFailed = relativeImbibedMassError > maxError;
    const bool saturationL1Failed = relativeSaturationL1Error > maxError;

    if (centerOfMassFailed || imbibedMassFailed || saturationL1Failed)
    {
        std::ostringstream message;
        message << "McWhorter semi-analytical check failed.";
        if (centerOfMassFailed)
            message << " Relative imbibed wetting-phase center-of-mass error "
                    << relativeCenterOfMassError << " exceeds threshold " << maxError << ".";
        if (imbibedMassFailed)
            message << " Relative imbibed wetting-phase mass error "
                    << relativeImbibedMassError << " exceeds threshold " << maxError << ".";
        if (saturationL1Failed)
            message << " Relative saturation L1 error "
                    << relativeSaturationL1Error << " exceeds threshold " << maxError << ".";
        DUNE_THROW(Dune::InvalidStateException, message.str());
    }

    if (leafGridView.comm().rank() == 0)
    {
        std::cout << "numeric center of imbibed wetting-phase mass: " << centerOfMassNumeric << std::endl;
        std::cout << "reference center of imbibed wetting-phase mass: " << centerOfMassReference << std::endl;
        std::cout << "numeric imbibed wetting-phase mass: " << imbibedMassNumeric << std::endl;
        std::cout << "reference imbibed wetting-phase mass: " << imbibedMassReference << std::endl;
        std::cout << "relative errors (mass, center of mass, saturation L1): "
                  << relativeImbibedMassError << ", "
                  << relativeCenterOfMassError << ", "
                  << relativeSaturationL1Error << std::endl;
    }

    if (leafGridView.comm().rank() == 0)
    {
        Parameters::print();
    }

    return 0;
}
