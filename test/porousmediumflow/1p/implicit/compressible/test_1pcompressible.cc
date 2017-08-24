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
 * \brief test for the one-phase CC model
 */
#include <config.h>

#include "problem.hh"

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/parameterparser.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonmethod.hh>

#include <dumux/implicit/cellcentered/assembler.hh>

#include <dumux/io/vtkoutputmodule.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = TTAG(IncompressibleTestProblem);

    // some aliases for better readability
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using ParameterTree = typename GET_PROP(TypeTag, ParameterTree);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);

    // TODO is needed for output. Should be in extra module
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    // for non-linear problems
    using NewtonController = typename GET_PROP_TYPE(TypeTag, NewtonController);

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // parse command line arguments
    ParameterParser::parseCommandLineArguments(argc, argv, ParameterTree::tree());

    // parse the input file into the parameter tree
    // check first if the user provided an input file through the command line, if not use the default
    const auto parameterFileName = ParameterTree::tree().hasKey("ParameterFile") ? GET_RUNTIME_PARAM(TypeTag, std::string, ParameterFile) : "";
    ParameterParser::parseInputFile(argc, argv, ParameterTree::tree(), parameterFileName);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    try { GridCreator::makeGrid(); }
    catch (...) {
        std::cout << "\n\t -> Creation of the grid failed! <- \n\n";
        throw;
    }
    GridCreator::loadBalance();

    // we compute on the leaf grid view
    const auto& leafGridView = GridCreator::grid().leafGridView();

    // create the finite volume grid geometry
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // the problem (initial and boundary conditions)
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    // the solution vector
    SolutionVector x(leafGridView.size(0));
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x, xOld);

    // read the time loop parameters
    auto tEnd = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeLoop, TEnd);
    auto dt = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeLoop, DtInitial);
    auto maxDivisions = GET_PARAM_FROM_GROUP(TypeTag, int, TimeLoop, MaxTimeStepDivisions);
    auto maxDt = GET_PARAM_FROM_GROUP(TypeTag, Scalar, TimeLoop, MaxTimeStepSize);

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (ParameterTree::tree().hasKey("Restart") || ParameterTree::tree().hasKey("TimeLoop.Restart"))
        restartTime = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeLoop, Restart);

    // write initial solution to disk
    // Dune::VTKSequenceWriter<GridView> vtkWriter(leafGridView, "test_1pcompressible", "", "");
    // vtkWriter.addCellData(x, "p");
    // vtkWriter.write(restartTime);

    // intialize the output module
    VtkOutputModule<TypeTag> vtkWriter(*problem, *fvGridGeometry, *gridVariables, x, problem->name());
    vtkWriter.addPrimaryVariable("pressure", Indices::pressureIdx);
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // make assembler
    auto assembler = std::make_shared<CCImplicitAssembler<TypeTag>>(problem, fvGridGeometry, gridVariables, timeLoop);

    // make linear solver
    auto linearSolver = std::make_shared<ILU0BiCGSTABBackend<TypeTag>>(*problem);

    // instantiate non-linear solver
    auto newtonController = std::make_shared<NewtonController>(leafGridView.comm(), timeLoop);
    NewtonMethod<TypeTag, NewtonController> nonLinearSolver(newtonController, assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // set solution
        assembler->setPreviousSolution(xOld);

        // try solving the non-linear system
        for (int i = 0; i < maxDivisions; ++i)
        {
            // linearize & solve
            auto converged = nonLinearSolver.solve(x);

            if (converged)
                break;

            if (!converged && i == maxDivisions-1)
                DUNE_THROW(Dune::MathError,
                           "Newton solver didn't converge after "
                           << maxDivisions
                           << " time-step divisions. dt="
                           << timeLoop->timeStepSize()
                           << ".\nThe solutions of the current and the previous time steps "
                           << "have been saved to restart files.");
        }

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();

        // write output
        vtkWriter.write(timeLoop->time());

        timeLoop->reportTimeStep();

        // set new dt
        timeLoop->setTimeStepSize(newtonController->suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    // print dumux end message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;

} // end main
