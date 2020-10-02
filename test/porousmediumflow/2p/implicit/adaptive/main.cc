// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup TwoPTests
 * \brief Test for the two-phase porous medium flow model
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_ug.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/adaptive/adapt.hh>
#include <dumux/adaptive/markelements.hh>
#include <dumux/adaptive/initializationindicator.hh>
#include <dumux/porousmediumflow/2p/griddatatransfer.hh>
#include <dumux/porousmediumflow/2p/gridadaptindicator.hh>

// Use the incompressible or point source problem for this adaptive test
#include <test/porousmediumflow/2p/implicit/incompressible/problem.hh>
#include "pointsourceproblem.hh"
#include "problem.hh"

// Type tags for the adaptive versions of the two-phase incompressible problem
namespace Dumux {
namespace Properties {
//! Type Tags for the adaptive tests
// Create new type tags
namespace TTag {
struct TwoPIncompressibleAdaptiveTpfa { using InheritsFrom = std::tuple<TwoPIncompressibleTpfa>; };
struct TwoPIncompressibleAdaptiveMpfa { using InheritsFrom = std::tuple<TwoPIncompressibleMpfa>; };
struct TwoPIncompressibleAdaptiveBox { using InheritsFrom = std::tuple<TwoPIncompressibleBox>; };
struct TwoPAdaptivePointSource { using InheritsFrom = std::tuple<TwoPIncompressibleAdaptiveTpfa>; };
} // end namespace TTag

//! Use non-conforming refinement in the cell-centered tests, conforming for box
#if HAVE_DUNE_ALUGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPIncompressibleAdaptiveTpfa> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPIncompressibleAdaptiveMpfa> { using type = Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>; };
#endif
#if HAVE_DUNE_UGGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::TwoPIncompressibleAdaptiveBox> { using type = Dune::UGGrid<2>; };
#endif
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPAdaptivePointSource> { using type = PointSourceTestProblem<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPIncompressibleAdaptiveTpfa> { using type = TwoPTestProblemAdaptive<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPIncompressibleAdaptiveMpfa> { using type = TwoPTestProblemAdaptive<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::TwoPIncompressibleAdaptiveBox> { using type = TwoPTestProblemAdaptive<TypeTag>; };
} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);
    problem->computePointSourceMap(); // enable point sources

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // instantiate indicator & data transfer, read parameters for indicator
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const Scalar refineTol = getParam<Scalar>("Adaptive.RefineTolerance");
    const Scalar coarsenTol = getParam<Scalar>("Adaptive.CoarsenTolerance");
    TwoPGridAdaptIndicator<TypeTag> indicator(gridGeometry);
    TwoPGridDataTransfer<TypeTag> dataTransfer(problem, gridGeometry, gridVariables, x);

    // Do initial refinement around sources/BCs
    GridAdaptInitializationIndicator<TypeTag> initIndicator(problem, gridGeometry, gridVariables);

    // refine up to the maximum level
    const auto maxLevel = getParam<std::size_t>("Adaptive.MaxLevel", 0);
    for (std::size_t i = 0; i < maxLevel; ++i)
    {
        initIndicator.calculate(x);

        bool wasAdapted = false;
        if (markElements(gridManager.grid(), initIndicator))
            wasAdapted = adapt(gridManager.grid(), dataTransfer);

        // update grid data after adaption
        if (wasAdapted)
        {
            xOld = x; //!< Overwrite the old solution with the new (resized & interpolated) one
            gridVariables->updateAfterGridAdaption(x); //!< Initialize the secondary variables to the new (and "new old") solution
            problem->computePointSourceMap(); //!< Update the point source map
        }
    }

    // Do refinement for the initial conditions using our indicator
    indicator.calculate(x, refineTol, coarsenTol);

    bool wasAdapted = false;
    if (markElements(gridManager.grid(), indicator))
        wasAdapted = adapt(gridManager.grid(), dataTransfer);

    // update grid data after adaption
    if (wasAdapted)
    {
        xOld = x; //!< Overwrite the old solution with the new (resized & interpolated) one
        gridVariables->updateAfterGridAdaption(x); //!< Initialize the secondary variables to the new (and "new old") solution
        problem->computePointSourceMap(); //!< Update the point source map
    }

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // refine/coarsen only after first time step
        if (timeLoop->time() > 0)
        {
            // compute refinement indicator
            indicator.calculate(x, refineTol, coarsenTol);

            // mark elements and maybe adapt grid
            wasAdapted = false;
            if (markElements(gridManager.grid(), indicator))
                wasAdapted = adapt(gridManager.grid(), dataTransfer);

            if (wasAdapted)
            {
                // Note that if we were using point sources, we would have to update the map here as well
                xOld = x; //!< Overwrite the old solution with the new (resized & interpolated) one
                assembler->setJacobianPattern(); //!< Tell the assembler to resize the matrix and set pattern
                assembler->setResidualSize(); //!< Tell the assembler to resize the residual
                gridVariables->updateAfterGridAdaption(x); //!< Initialize the secondary variables to the new (and "new old") solution
                problem->computePointSourceMap(); //!< Update the point source map
            }
        }

        // solve the non-linear system with time step control
        nonLinearSolver.solve(x, *timeLoop);

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

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
