// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVETests
 * \brief Integration test for the immiscible two-phase vertical-equilibrium model in two and three dimensions.
 *
 * The test solves a heterogeneous injection problem on the vertically
 * integrated coarse grid and reconstructs fine-level quantities in each
 * vertical column. See the documentation of the TwoPVE model for the
 * assumptions and restrictions of the current implementation.
 */

#include <config.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <memory>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/porousmediumflow/2pve/finelevelview.hh>

#include "properties.hh"
#include "massbalance.hh"

#ifndef DIFFMETHOD
#define DIFFMETHOD DiffMethod::numeric
#endif

namespace Dumux::VETest {

/*!
 * \brief Initializes the horizontal coarse grid by projecting the fine-grid bounds.
 *
 * \param gridManagerCoarse gridManager belonging to the coarse grid
 */
template<class GridManager, class TypeTag>
void initializeCoarseGrid(GridManager& gridManagerCoarse)
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    constexpr int dim = GridView::dimension;

    using FineGrid = GetPropType<TypeTag, Properties::FineGrid>;
    using FinePosition = Dune::FieldVector<double, FineGrid::dimensionworld>;
    const auto lowerLeftFine = getParam<FinePosition>("Grid.LowerLeft");
    const auto upperRightFine = getParam<FinePosition>("Grid.UpperRight");
    const auto cellsFine = getParam<std::array<int, FineGrid::dimension>>("Grid.Cells");
    GlobalPosition lowerLeftCoarse, upperRightCoarse;
    std::array<int, dim> cellsCoarse;
    for (int axis = 0; axis < dim; ++axis)
    {
        lowerLeftCoarse[axis] = lowerLeftFine[axis];
        upperRightCoarse[axis] = upperRightFine[axis];
        cellsCoarse[axis] = cellsFine[axis];
    }
    gridManagerCoarse.init(lowerLeftCoarse, upperRightCoarse, cellsCoarse);
}

/*!
 * \brief Concatenates the coarse-level output name
 *
 * \param problemCoarse coarse-level problem
 */
template<typename TypeTag>
std::string outputNameCoarse(const GetPropType<TypeTag, Properties::Problem>& problemCoarse)
{
    using GridView = typename GetPropType<TypeTag, Properties::FineGridGeometry>::GridView;
    constexpr int dim = GridView::dimension;
    using CellArray = std::array<unsigned int, GridView::dimensionworld>;
    const auto numberCellsFine = getParam<CellArray>("Grid.Cells");
    const std::size_t numCellsDim1 = numberCellsFine[0];
    std::size_t numCellsDim2 = 0;
    if constexpr (dim == 3)
        numCellsDim2 = numberCellsFine[1];
    const std::size_t numCellsDim3 = numberCellsFine[dim-1];

    return problemCoarse.name() + std::to_string(numCellsDim1) + "x" + std::to_string(numCellsDim2) + "x" + std::to_string(numCellsDim3);
}

} // end namespace Dumux::VETest

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TwoPVEImmiscibleTpfa;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // create a grid
    using GridManagerType = GridManager<GetPropType<TypeTag, Properties::Grid>>;
    GridManager<GetPropType<TypeTag, Properties::FineGrid>> gridManagerFine;
    gridManagerFine.init("");
    GridManagerType gridManagerCoarse;
    Dumux::VETest::initializeCoarseGrid<GridManagerType, TypeTag>(gridManagerCoarse);

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridViewCoarse = gridManagerCoarse.grid().leafGridView();
    const auto& leafGridViewFine = gridManagerFine.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometryCoarse = std::make_shared<GridGeometry>(leafGridViewCoarse);
    using FineGridGeometry = GetPropType<TypeTag, Properties::FineGridGeometry>;
    auto gridGeometryFine = std::make_shared<FineGridGeometry>(leafGridViewFine);

    // create the fine-level problem and the fine-level view of the VE scheme
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FineProblem = typename Problem::FineProblem;
    using FineLevelView = typename Problem::FineLevelView;
    auto fineProblem = std::make_shared<FineProblem>(gridGeometryFine);
    auto fineLevelView = std::make_shared<FineLevelView>(gridGeometryFine, gridGeometryCoarse, fineProblem->spatialParamsPtr());

    // the problem (initial and boundary conditions)
    auto problemCoarse = std::make_shared<Problem>(gridGeometryCoarse, fineLevelView, fineProblem);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    const bool usePrescribedTimeSteps = getParam<bool>("TimeLoop.UsePrescribedTimeSteps", false);
    // Keep the nominal ramp independent of Newton iterations and checkpoint clipping.
    Scalar prescribedDt = std::min(dt, maxDt);
    std::size_t retriedTimeSteps = 0;

    // instantiate time loop
    auto timeLoopCoarse = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoopCoarse->setMaxTimeStepSize(maxDt);
    timeLoopCoarse->setPeriodicCheckPoint(tEnd/10.0);
    problemCoarse->setTimeLoop(timeLoopCoarse);
    if (problemCoarse->injectionEndTime() < tEnd)
        timeLoopCoarse->setCheckPoint(problemCoarse->injectionEndTime());

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    auto xCoarse = std::make_shared<SolutionVector>(gridGeometryCoarse->numDofs());

    // initialize solution
    problemCoarse->applyInitialSolution(*xCoarse);
    fineLevelView->updateSol(*problemCoarse, *xCoarse);

    auto xOldCoarse = *xCoarse;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariablesCoarse = std::make_shared<GridVariables>(problemCoarse, gridGeometryCoarse);
    gridVariablesCoarse->init(*xCoarse);

    // coarse-level vtk output
    const std::string coarseOutputName = Dumux::VETest::outputNameCoarse<TypeTag>(*problemCoarse);
    VtkOutputModule<GridVariables, SolutionVector> vtkWriterCoarse(*gridVariablesCoarse, *xCoarse, coarseOutputName);
    using IOFieldsVECoarse = GetPropType<TypeTag, Properties::IOFields>;
    using VelocityOutputVE = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriterCoarse.addVelocityOutput(std::make_shared<VelocityOutputVE>(*gridVariablesCoarse));
    vtkWriterCoarse.addVolumeVariable([](const auto& v){return v.gasPlumeDist();}, "zp");
    IOFieldsVECoarse::initOutputModule(vtkWriterCoarse);
    vtkWriterCoarse.addVolumeVariable([](const auto& v){ return v.porosityIntegral(); }, "porosityIntegral");
    vtkWriterCoarse.addVolumeVariable([](const auto& v){ return v.permeabilityIntegral(); }, "permeabilityIntegral");
    vtkWriterCoarse.write(0.0);

    // fine-level vtk output
    using GridView = typename FineGridGeometry::GridView;
    const std::string fineOutputName = "fine_" + Dumux::VETest::outputNameCoarse<TypeTag>(*problemCoarse);
    Dune::VTKSequenceWriter<GridView> vtkWriterFineLevel(gridGeometryFine->gridView(), fineOutputName, ".", "");
    fineLevelView->fields().registerFields(vtkWriterFineLevel);
    vtkWriterFineLevel.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DIFFMETHOD>;
    auto assemblerCoarse = std::make_shared<Assembler>(problemCoarse, gridGeometryCoarse, gridVariablesCoarse, timeLoopCoarse, xOldCoarse);

    // the linear solver
    using LinearSolver = ILUBiCGSTABIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolverCoarse = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolverCoarse(assemblerCoarse, linearSolverCoarse);

    const auto massBalanceTolerance = getParam<Scalar>("MassBalance.RelativeTolerance");

    // time loop (solution is conducted on coarse level of VE scheme)
    timeLoopCoarse->start(); do
    {
        // solve the non-linear system with time step control
        const auto requestedDt = timeLoopCoarse->timeStepSize();
        nonLinearSolverCoarse.solve(*xCoarse, *timeLoopCoarse);
        if (timeLoopCoarse->timeStepSize() < requestedDt)
        {
            ++retriedTimeSteps;
            if (mpiHelper.rank() == 0)
                std::cout << "Newton retry reduced time step from " << requestedDt
                          << " to " << timeLoopCoarse->timeStepSize() << " seconds\n";
        }

        // the cached coarse-level volume variables depend on the column history updated with the fine-level solution
        fineLevelView->updateSol(*problemCoarse, *xCoarse);
        gridVariablesCoarse->update(*xCoarse);

        // make the new solution the old solution
        xOldCoarse = *xCoarse;
        gridVariablesCoarse->advanceTimeStep();

        // advance the time loop to the next step
        timeLoopCoarse->advanceTimeStep();

        // compute and print mass balance for gas phase
        const auto massBalance = Dumux::VETest::computeMassBalance<TypeTag>(*gridGeometryCoarse, *gridGeometryFine, *xCoarse, *problemCoarse, *timeLoopCoarse);
        Dumux::VETest::printMassBalance(massBalance);
        Dumux::VETest::checkMassBalance(massBalance, massBalanceTolerance);

        // write vtk files for coarse and fine level
        if (timeLoopCoarse->isCheckPoint() || timeLoopCoarse->finished())
        {
            vtkWriterCoarse.write(timeLoopCoarse->time());
            vtkWriterFineLevel.write(timeLoopCoarse->time());
        }

        // report statistics of this time step
        timeLoopCoarse->reportTimeStep();

        if (usePrescribedTimeSteps)
        {
            prescribedDt = std::min(2.0*prescribedDt, maxDt);
            timeLoopCoarse->setTimeStepSize(prescribedDt);
        }
        else
            timeLoopCoarse->setTimeStepSize(nonLinearSolverCoarse.suggestTimeStepSize(timeLoopCoarse->timeStepSize()));
    } while (!timeLoopCoarse->finished());

    const auto finalMassBalance = Dumux::VETest::computeMassBalance<TypeTag>(*gridGeometryCoarse, *gridGeometryFine, *xCoarse, *problemCoarse, *timeLoopCoarse);
    const auto minimumTrappedFraction = getParam<Scalar>("MassBalance.MinTrappedGasFraction", 0.0);
    const Scalar trappedFraction = finalMassBalance.trappedGasMass/finalMassBalance.expectedInjectedMass;
    if (mpiHelper.rank() == 0)
        std::cout << "Final residually trapped gas mass: " << finalMassBalance.trappedGasMass
                  << (GridView::dimension == 2 ? " kg per metre out-of-plane" : " kg")
                  << " (" << 100*trappedFraction << "% of injected gas)\n";
    if (trappedFraction < minimumTrappedFraction)
        DUNE_THROW(Dune::Exception, "Trapped gas fraction " << trappedFraction
                   << " is below the required minimum " << minimumTrappedFraction);

    // output some Newton statistics
    if (usePrescribedTimeSteps && mpiHelper.rank() == 0)
        std::cout << "Prescribed schedule: " << retriedTimeSteps << " time steps required Newton retries\n";
    nonLinearSolverCoarse.report();

    timeLoopCoarse->finalize(leafGridViewCoarse.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
