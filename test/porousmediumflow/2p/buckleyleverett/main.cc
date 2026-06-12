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


    using ScalarVector = typename BuckleyLeverettAnalyticSolution<TypeTag>::ScalarVector;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    constexpr auto saturationIdx = ModelTraits::Indices::saturationIdx;
    using BarycenterVector = Dune::FieldVector<Scalar, 2>;

    // compute relative error in wetting-phase mass barycenters
    // assumptions: solutions are (pseudo)1D, densities and porosity are constant, TPFA discretization
    auto computeBarycentersMassWX = [&problem, &gridGeometry](const SolutionVector& numericSolution,
                                                              const ScalarVector& analyticSolution) -> BarycenterVector
    {
        Scalar leftXBoundary = getParam<GlobalPosition>("Grid.LowerLeft")[0];
        FluidState fluidState;
        Scalar referencePressure = getParam<Scalar>("Problem.ReferencePressure");
        fluidState.setTemperature(problem->spatialParams().temperatureAtPos(GlobalPosition{}));
        fluidState.setPressure(FluidSystem::phase0Idx, referencePressure);
        fluidState.setPressure(FluidSystem::phase1Idx, referencePressure);

        const auto densityW = FluidSystem::density(fluidState, FluidSystem::phase0Idx);
        const auto porosity = problem->spatialParams().porosityAtPos(GlobalPosition{});

        Scalar barycenterMassWNumeric = 0.0;
        Scalar barycenterMassWAnalytic = 0.0;
        Scalar totalMassWNumeric = 0.0;
        Scalar totalMassWAnalytic = 0.0;
        BarycenterVector barycenters(0.0);

        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto globalPos = element.geometry().center();
            const auto eIdx = gridGeometry->elementMapper().index(element);
            const auto volume = element.geometry().volume();

            Scalar satWNumeric = 1.0-numericSolution[eIdx][saturationIdx];
            Scalar satWAnalytic = analyticSolution[eIdx];

            Scalar localMassWNumeric = densityW * porosity * volume * satWNumeric;
            Scalar localMassWAnalytic = densityW * porosity * volume * satWAnalytic;

            totalMassWNumeric += localMassWNumeric;
            totalMassWAnalytic += localMassWAnalytic;

            barycenterMassWNumeric += globalPos[0] * localMassWNumeric;
            barycenterMassWAnalytic += globalPos[0] * localMassWAnalytic;
        }

        // avoid dividing by zero
        if(totalMassWNumeric < 1e-16)
            barycenterMassWNumeric = leftXBoundary;
        else
            barycenterMassWNumeric /= totalMassWNumeric;

        if(totalMassWAnalytic < 1e-16)
            barycenterMassWAnalytic = leftXBoundary;
        else
            barycenterMassWAnalytic /= totalMassWAnalytic;

        barycenters[0] = barycenterMassWNumeric;
        barycenters[1] = barycenterMassWAnalytic;

        return barycenters;
    };

    auto finalBarycentersX = computeBarycentersMassWX(x, analyticSolution.values());

    const auto maxBarycenterMassWError = getParam<Scalar>("Problem.MaxBarycenterMassWError");
    auto barycenterMassWError = 1.0;
    // avoid dividing by zero
    if(std::abs(finalBarycentersX[1])<1e-16)
        barycenterMassWError = finalBarycentersX[1];
    else
        barycenterMassWError = std::abs(finalBarycentersX[1] - finalBarycentersX[0])/finalBarycentersX[1];
    if (barycenterMassWError > maxBarycenterMassWError)
        DUNE_THROW(Dune::InvalidStateException, "Relative wetting-phase mass barycenter error " << barycenterMassWError
                   << " exceeds the threshold " << maxBarycenterMassWError);

    if (leafGridView.comm().rank() == 0)
    {
        std::cout << "Relative wetting-phase mass barycenter error against analytical Buckley-Leverett solution: "
                  << barycenterMassWError << std::endl;
        Parameters::print();
    }

    return 0;
}
