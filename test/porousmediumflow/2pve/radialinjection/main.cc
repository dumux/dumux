// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
 * \brief Radial injection into a homogeneous confined aquifer compared to the similarity solution of Nordbotten and Celia (2006).
 */

#include <config.h>

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "analyticsolution.hh"
#include "properties.hh"

namespace Dumux {

/*!
 * \brief The mean deviation of the interface height from the similarity solution over the plume extent, relative to the aquifer height
 */
template<class GridGeometry, class Scalar>
Scalar interfaceError(const GridGeometry& gridGeometry,
                      const std::vector<Scalar>& interfaceHeight,
                      const std::vector<Scalar>& exactInterfaceHeight,
                      Scalar aquiferHeight,
                      Scalar plumeExtent)
{
    using std::abs;
    Scalar error = 0.0;
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto elementIdx = gridGeometry.elementMapper().index(element);
        const auto& geometry = element.geometry();
        const Scalar radialExtent = geometry.corner(1)[0] - geometry.corner(0)[0];
        error += abs(interfaceHeight[elementIdx] - exactInterfaceHeight[elementIdx])*radialExtent;
    }
    return error/(aquiferHeight*plumeExtent);
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TwoPVERadialInjection;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);

    // the fine grid resolves the vertical direction, the coarse grid consists of a single layer of columns
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GlobalPosition = Dune::FieldVector<double, Grid::dimensionworld>;
    GridManager<Grid> gridManagerFine;
    gridManagerFine.init();
    GridManager<Grid> gridManagerCoarse;
    auto coarseCells = getParam<std::array<int, Grid::dimension>>("Grid.Cells");
    coarseCells[Grid::dimension-1] = 1;
    gridManagerCoarse.init(getParam<GlobalPosition>("Grid.LowerLeft"), getParam<GlobalPosition>("Grid.UpperRight"), coarseCells);

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometryCoarse = std::make_shared<GridGeometry>(gridManagerCoarse.grid().leafGridView());
    auto gridGeometryFine = std::make_shared<GridGeometry>(gridManagerFine.grid().leafGridView());

    // the fine-level view and the coarse-level spatial parameters upscaled from the fine level
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;
    using FineLevelView = typename Problem::FineLevelView;
    auto spatialParamsFine = std::make_shared<typename SpatialParams::SpatialParamsFine>(gridGeometryFine);
    auto fineLevelView = std::make_shared<FineLevelView>(gridGeometryFine, gridGeometryCoarse, spatialParamsFine);
    auto spatialParams = std::make_shared<SpatialParams>(gridGeometryCoarse, fineLevelView->columnMap(), spatialParamsFine, fineLevelView->fineCellHeight());
    auto problem = std::make_shared<Problem>(gridGeometryCoarse, spatialParams, fineLevelView);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometryCoarse->numDofs());
    problem->applyInitialSolution(x);
    fineLevelView->updateSol(*problem, x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometryCoarse);
    gridVariables->init(x);

    // the similarity solution for the mobility ratio of the sharp-interface limit, in which both phases flow with unit relative permeability
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    GetPropType<TypeTag, Properties::FluidState> fluidState;
    fluidState.setTemperature(spatialParams->temperatureAtPos(gridGeometryCoarse->bBoxMin()));
    fluidState.setPressure(FluidSystem::phase0Idx, getParam<Scalar>("Problem.InitialPressure"));
    fluidState.setPressure(FluidSystem::phase1Idx, getParam<Scalar>("Problem.InitialPressure"));
    const Scalar viscosityResident = FluidSystem::viscosity(fluidState, FluidSystem::phase0Idx);
    const Scalar viscosityInjected = FluidSystem::viscosity(fluidState, FluidSystem::phase1Idx);
    const Scalar densityDifference = FluidSystem::density(fluidState, FluidSystem::phase0Idx) - FluidSystem::density(fluidState, FluidSystem::phase1Idx);
    const Scalar aquiferHeight = gridGeometryCoarse->bBoxMax()[1] - gridGeometryCoarse->bBoxMin()[1];
    const Scalar permeability = getParam<Scalar>("SpatialParams.Permeability");
    const TwoPVERadialInjectionSimilaritySolution<Scalar> similaritySolution(viscosityResident/viscosityInjected,
                                                                            aquiferHeight,
                                                                            getParam<Scalar>("SpatialParams.Porosity"),
                                                                            getParam<Scalar>("SpatialParams.Swr"),
                                                                            problem->injectionRate());
    const Scalar gravityNumber = 2.0*M_PI*densityDifference*spatialParams->gravity(gridGeometryCoarse->bBoxMin()).two_norm()
                                 *permeability/viscosityResident*aquiferHeight*aquiferHeight/problem->injectionRate();
    std::cout << "Mobility ratio: " << viscosityResident/viscosityInjected << ", gravity number: " << gravityNumber << std::endl;

    // the effective interface height corresponds to a plume that contains the injected fluid at the saturation 1 - Swr
    const Scalar residualSaturation = getParam<Scalar>("SpatialParams.Swr");
    std::vector<Scalar> exactInterfaceHeight(gridGeometryCoarse->gridView().size(0));
    std::vector<Scalar> interfaceHeight(gridGeometryCoarse->gridView().size(0));
    const auto updateInterfaceHeights = [&](Scalar time)
    {
        for (const auto& element : elements(gridGeometryCoarse->gridView()))
        {
            const auto elementIdx = gridGeometryCoarse->elementMapper().index(element);
            const Scalar radius = element.geometry().center()[0];
            exactInterfaceHeight[elementIdx] = time > 0.0 ? similaritySolution.interfaceHeight(radius, time) : aquiferHeight;
            const Scalar saturationInjected = gridVariables->curGridVolVars().volVars(elementIdx).saturation(FluidSystem::phase1Idx);
            interfaceHeight[elementIdx] = aquiferHeight*(1.0 - saturationInjected/(1.0 - residualSaturation));
        }
    };

    using VtkOutputFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    VtkOutputFields::initOutputModule(vtkWriter);
    vtkWriter.addVolumeVariable([](const auto& v){ return v.gasPlumeDist(); }, "zp");
    vtkWriter.addField(interfaceHeight, "interfaceHeight");
    vtkWriter.addField(exactInterfaceHeight, "interfaceHeightExact");
    updateInterfaceHeights(0.0);
    vtkWriter.write(0.0);

    using GridView = typename GridGeometry::GridView;
    Dune::VTKSequenceWriter<GridView> vtkWriterFine(gridGeometryFine->gridView(), "fine_" + problem->name(), ".", "");
    fineLevelView->fields().registerFields(vtkWriterFine);
    vtkWriterFine.write(0.0);

    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, getParam<Scalar>("TimeLoop.DtInitial"), tEnd);
    timeLoop->setMaxTimeStepSize(getParam<Scalar>("TimeLoop.MaxTimeStepSize"));
    timeLoop->setPeriodicCheckPoint(tEnd/getParam<int>("TimeLoop.NumOutputs"));

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometryCoarse, gridVariables, timeLoop, xOld);
    using LinearSolver = ILUBiCGSTABIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    Scalar error = 0.0;
    timeLoop->start(); do
    {
        nonLinearSolver.solve(x, *timeLoop);

        // the cached coarse-level volume variables depend on the column history updated with the fine-level solution
        fineLevelView->updateSol(*problem, x);
        gridVariables->update(x);

        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        if (timeLoop->isCheckPoint() || timeLoop->finished())
        {
            updateInterfaceHeights(timeLoop->time());
            error = interfaceError(*gridGeometryCoarse, interfaceHeight, exactInterfaceHeight, aquiferHeight, similaritySolution.plumeExtent(timeLoop->time()));
            std::cout << "Relative interface error at t = " << timeLoop->time() << " s: " << error << std::endl;
            vtkWriter.write(timeLoop->time());
            vtkWriterFine.write(timeLoop->time());
        }

        timeLoop->reportTimeStep();
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
    } while (!timeLoop->finished());

    nonLinearSolver.report();
    timeLoop->finalize(gridGeometryCoarse->gridView().comm());

    const Scalar maxError = getParam<Scalar>("Benchmark.MaxRelativeInterfaceError");
    if (error > maxError)
        DUNE_THROW(Dune::Exception, "Relative interface error " << error << " exceeds " << maxError);

    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
