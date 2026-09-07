// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup VETest
 * \brief Integration test for the immiscible two-phase vertical-equilibrium model in two and three dimensions.
 *
 * The test solves a heterogeneous injection problem on the vertically
 * integrated coarse grid and reconstructs fine-level quantities in each
 * vertical column. See the documentation of the TwoPVE model for the
 * assumptions and restrictions of the current implementation.
 */

#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/porousmediumflow/2pve/finelevel_view.hh>

#include "properties.hh"
#include "massbalance.hh"
#include "problem_fine.hh"

#ifndef DIFFMETHOD
#define DIFFMETHOD DiffMethod::numeric
#endif

namespace Dumux::VETest {

/*!
 * \brief Initializes the coarse grid. Sets the number of vertical cells to 1.
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

    const GlobalPosition lowerLeftCoarse = getParam<GlobalPosition>("Grid.LowerLeft");
    const GlobalPosition upperRightCoarse = getParam<GlobalPosition>("Grid.UpperRight");
    std::array<int, dim> cellsCoarse = getParam<std::array<int, dim>>("Grid.Cells");
    cellsCoarse[dim-1] = 1;
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
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    constexpr int dim = GridView::dimension;
    using CellArray = std::array<unsigned int, GridView::dimensionworld>;
    const auto numberCellsFine = getParam<CellArray>("Grid.Cells");
    size_t numCellsDim1 = numberCellsFine[0];
    size_t numCellsDim2 = 0;
    if constexpr(dim==3)
        numCellsDim2 = numberCellsFine[1];
    size_t numCellsDim3 = numberCellsFine[dim-1];

    return problemCoarse.name() + std::to_string(numCellsDim1) + "x" + std::to_string(numCellsDim2) + "x" + std::to_string(numCellsDim3);
}

} // end namespace Dumux::VETest


int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TwoPVEImmiscibleTpfa;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // create a grid
    using GridManagerType = GridManager<GetPropType<TypeTag, Properties::Grid>>;
    GridManagerType gridManagerFine;
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
    auto gridGeometryFine = std::make_shared<GridGeometry>(leafGridViewFine);

    // create fine-level view of VE scheme
    using FineProblem = TwoPVEFineProblem<TypeTag>;
    using FineLevelView = TwoPVEFineLevelView<TypeTag, FineProblem>;
    auto fineLevelView = std::make_shared<FineLevelView>(gridGeometryFine, gridGeometryCoarse);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problemCoarse = std::make_shared<Problem>(gridGeometryCoarse, fineLevelView);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // instantiate time loop
    auto timeLoopCoarse = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoopCoarse->setMaxTimeStepSize(maxDt);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    auto xCoarse = std::make_shared<SolutionVector>(gridGeometryCoarse->numDofs());

    // intialize solution
    problemCoarse->applyInitialSolution(*xCoarse);
    fineLevelView->updateSol(*problemCoarse, *xCoarse);

    auto xOldCoarse = *xCoarse;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariablesCoarse = std::make_shared<GridVariables>(problemCoarse, gridGeometryCoarse);
    gridVariablesCoarse->init(*xCoarse);

    // coarse-level vtk output
    std::string VECoarseOutputName = Dumux::VETest::outputNameCoarse<TypeTag>(*problemCoarse);
    VtkOutputModule<GridVariables, SolutionVector> vtkWriterCoarse(*gridVariablesCoarse, *xCoarse, VECoarseOutputName);
    using IOFieldsVECoarse = GetPropType<TypeTag, Properties::IOFields>;
    using VelocityOutputVE = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriterCoarse.addVelocityOutput(std::make_shared<VelocityOutputVE>(*gridVariablesCoarse));
    vtkWriterCoarse.addVolumeVariable([](const auto& v){return v.permeability();}, "permeability");
    vtkWriterCoarse.addVolumeVariable([](const auto& v){return v.gasPlumedist();}, "zp");
    IOFieldsVECoarse::initOutputModule(vtkWriterCoarse);
    vtkWriterCoarse.write(0.0);

    // fine-level vtk output
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    std::string VEFineOutputName = "fine_" + Dumux::VETest::outputNameCoarse<TypeTag>(*problemCoarse);
    Dune::VTKSequenceWriter<GridView> vtkWriterFineLevel(gridGeometryFine->gridView(), VEFineOutputName, ".", "");
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

    // time loop (solution is conducted on coarse level of VE scheme)
    timeLoopCoarse->start(); do
    {
        // solve the non-linear system with time step control
        nonLinearSolverCoarse.solve(*xCoarse, *timeLoopCoarse);

        // const auto coarseBefore = *xCoarse;
        fineLevelView->updateSol(*problemCoarse, *xCoarse);

        // make the new solution the old solution
        xOldCoarse = *xCoarse;
        gridVariablesCoarse->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoopCoarse->advanceTimeStep();

        // compute and print mass balance for gas phase
        const auto massBalance = Dumux::VETest::computeMassBalance<TypeTag>(*gridGeometryCoarse, *gridGeometryFine, *xCoarse, *problemCoarse, *timeLoopCoarse);
        Dumux::VETest::printMassBalance(massBalance);

        // write vtk files for coarse and fine level
        vtkWriterCoarse.write(timeLoopCoarse->time());
        vtkWriterFineLevel.write(timeLoopCoarse->time());

        // report statistics of this time step
        timeLoopCoarse->reportTimeStep();

        // set new dt as suggested by the Newton solver
        timeLoopCoarse->setTimeStepSize(nonLinearSolverCoarse.suggestTimeStepSize(timeLoopCoarse->timeStepSize()));
    } while (!timeLoopCoarse->finished());

    // output some Newton statistics
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
