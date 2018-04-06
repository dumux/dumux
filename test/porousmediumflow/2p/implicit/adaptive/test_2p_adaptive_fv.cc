// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief test for the two-phase porousmedium flow model
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/methods.hh>

#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/adaptive/adapt.hh>
#include <dumux/adaptive/markelements.hh>
#include <dumux/adaptive/initializationindicator.hh>
#include <dumux/porousmediumflow/2p/griddatatransfer.hh>
#include <dumux/porousmediumflow/2p/gridadaptindicator.hh>

//! Use the incompressible or point source problem for this adaptive test
#include <test/porousmediumflow/2p/implicit/incompressible/problem.hh>
#include "pointsourceproblem.hh"
#include "problem.hh"

//! type tags for the adaptive versions of the two-phase incompressible problem
namespace Dumux {
namespace Properties {
//! Type Tags for the adaptive tests
NEW_TYPE_TAG(TwoPIncompressibleAdaptiveTpfa, INHERITS_FROM(TwoPIncompressibleTpfa));
NEW_TYPE_TAG(TwoPIncompressibleAdaptiveMpfa, INHERITS_FROM(TwoPIncompressibleMpfa));
NEW_TYPE_TAG(TwoPIncompressibleAdaptiveBox, INHERITS_FROM(TwoPIncompressibleBox));
NEW_TYPE_TAG(TwoPAdaptivePointSource, INHERITS_FROM(TwoPIncompressibleAdaptiveTpfa));

//! Use non-conforming refinement in the cell-centered tests, conforming for box
#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(TwoPIncompressibleAdaptiveTpfa, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
SET_TYPE_PROP(TwoPIncompressibleAdaptiveMpfa, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
#endif
#if HAVE_DUNE_UGGRID
SET_TYPE_PROP(TwoPIncompressibleAdaptiveBox, Grid, Dune::UGGrid<2>);
#endif
SET_TYPE_PROP(TwoPAdaptivePointSource, Problem, PointSourceTestProblem<TypeTag>);
SET_TYPE_PROP(TwoPIncompressibleAdaptiveTpfa, Problem, TwoPTestProblemAdaptive<TypeTag>);
SET_TYPE_PROP(TwoPIncompressibleAdaptiveMpfa, Problem, TwoPTestProblemAdaptive<TypeTag>);
SET_TYPE_PROP(TwoPIncompressibleAdaptiveBox, Problem, TwoPTestProblemAdaptive<TypeTag>);
} // end namespace Properties
} // end namespace Dumux

int main(int argc, char** argv) try
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = TTAG(TYPETAG);

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    GridCreator::makeGrid();
    GridCreator::loadBalance();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = GridCreator::grid().leafGridView();

    // create the finite volume grid geometry
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    auto problem = std::make_shared<Problem>(fvGridGeometry);
    problem->computePointSourceMap(); // enable point sources

    // the solution vector
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    SolutionVector x(fvGridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x, xOld);

    // instantiate indicator & data transfer, read parameters for indicator
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    const Scalar refineTol = getParam<Scalar>("Adaptive.RefineTolerance");
    const Scalar coarsenTol = getParam<Scalar>("Adaptive.CoarsenTolerance");
    TwoPGridAdaptIndicator<TypeTag> indicator(fvGridGeometry);
    TwoPGridDataTransfer<TypeTag> dataTransfer(problem, fvGridGeometry, gridVariables, x);

    // Do initial refinement around sources/BCs
    GridAdaptInitializationIndicator<TypeTag> initIndicator(problem, fvGridGeometry, gridVariables);

    // refine up to the maximum level
    const auto maxLevel = getParam<std::size_t>("Adaptive.MaxLevel", 0);
    for (std::size_t i = 0; i < maxLevel; ++i)
    {
        initIndicator.calculate(x);

        bool wasAdapted = false;
        if (markElements(GridCreator::grid(), initIndicator))
            wasAdapted = adapt(GridCreator::grid(), dataTransfer);

        // update grid data after adaption
        if (wasAdapted)
        {
            xOld = x;                         //!< Overwrite the old solution with the new (resized & interpolated) one
            gridVariables->init(x, xOld);     //!< Initialize the secondary variables to the new (and "new old") solution
            problem->computePointSourceMap(); //!< Update the point source map
        }
    }

    // Do refinement for the initial conditions using our indicator
    indicator.calculate(x, refineTol, coarsenTol);

    bool wasAdapted = false;
    if (markElements(GridCreator::grid(), indicator))
        wasAdapted = adapt(GridCreator::grid(), dataTransfer);

    // update grid data after adaption
    if (wasAdapted)
    {
        xOld = x;                         //!< Overwrite the old solution with the new (resized & interpolated) one
        gridVariables->init(x, xOld);     //!< Initialize the secondary variables to the new (and "new old") solution
        problem->computePointSourceMap(); //!< Update the point source map
    }

    // get some time loop parameters
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (haveParam("Restart") || haveParam("TimeLoop.Restart"))
        restartTime = getParam<Scalar>("TimeLoop.Restart");

    // intialize the vtk output module
    using VtkOutputFields = typename GET_PROP_TYPE(TypeTag, VtkOutputFields);
    VtkOutputModule<TypeTag> vtkWriter(*problem, *fvGridGeometry, *gridVariables, x, problem->name());
    VtkOutputFields::init(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop);

    // the linear solver
    using LinearSolver = AMGBackend<TypeTag>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, fvGridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // refine/coarsen only after first time step
        if (timeLoop->time() > restartTime)
        {
            // compute refinement indicator
            indicator.calculate(x, refineTol, coarsenTol);

            // mark elements and maybe adapt grid
            bool wasAdapted = false;
            if (markElements(GridCreator::grid(), indicator))
                wasAdapted = adapt(GridCreator::grid(), dataTransfer);

            if (wasAdapted)
            {
                // Note that if we were using point sources, we would have to update the map here as well
                xOld = x;                         //!< Overwrite the old solution with the new (resized & interpolated) one
                assembler->setJacobianPattern();  //!< Tell the assembler to resize the matrix and set pattern
                assembler->setResidualSize();     //!< Tell the assembler to resize the residual
                gridVariables->init(x, xOld);     //!< Initialize the secondary variables to the new (and "new old") solution
                problem->computePointSourceMap(); //!< Update the point source map
            }
        }

        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

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
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
