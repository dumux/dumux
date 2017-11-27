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
 * \brief test for the tracer CC model
 */
#include <config.h>

#include "tracertestproblem.hh"

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

#include <dumux/assembly/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>

int main(int argc, char** argv) try
{
    using namespace Dumux;

    //! define the type tag for this problem
    using TypeTag = TTAG(TYPETAG);

    //! initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    //! print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    //! parse command line arguments
    using ParameterTree = typename GET_PROP(TypeTag, ParameterTree);
    ParameterParser::parseCommandLineArguments(argc, argv, ParameterTree::tree());

    //! parse the input file into the parameter tree
    //! check first if the user provided an input file through the command line, if not use the default
    const auto parameterFileName = ParameterTree::tree().hasKey("ParameterFile") ?
        GET_RUNTIME_PARAM(TypeTag, std::string, ParameterFile) : "";
    ParameterParser::parseInputFile(argc, argv, ParameterTree::tree(), parameterFileName);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    try { GridCreator::makeGrid(); }
    catch (...) {
        std::cout << "\n\t -> Creation of the grid failed! <- \n\n";
        throw;
    }
    GridCreator::loadBalance();

    ////////////////////////////////////////////////////////////
    // setup instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    //! we compute on the leaf grid view
    const auto& leafGridView = GridCreator::grid().leafGridView();

    //! create the finite volume grid geometry
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    //! the problem (initial and boundary conditions)
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    //! the solution vector
    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    static constexpr int dofCodim = isBox ? GridView::dimension : 0;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    SolutionVector x(leafGridView.size(dofCodim));
    problem->applyInitialSolution(x);
    auto xOld = x;

    //! the grid variables
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x, xOld);

    //! get some time loop parameters
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    auto tEnd = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeLoop, TEnd);
    auto dt = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeLoop, DtInitial);
    auto maxDt = GET_PARAM_FROM_GROUP(TypeTag, Scalar, TimeLoop, MaxTimeStepSize);

    //! check if we are about to restart a previously interrupted simulation
    Scalar restartTime = 0;
    if (ParameterTree::tree().hasKey("Restart") || ParameterTree::tree().hasKey("TimeLoop.Restart"))
        restartTime = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeLoop, Restart);

    //! instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    //! the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic, IMPLICIT>;
    auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<SolutionVector>();
    assembler->setLinearSystem(A, r);

    //! the linear solver
    using LinearSolver = UMFPackBackend<TypeTag>;
    auto linearSolver = std::make_shared<LinearSolver>(*problem);

    //! intialize the vtk output module
    VtkOutputModule<TypeTag> vtkWriter(*problem, *fvGridGeometry, *gridVariables, x, problem->name());
    using VtkOutputFields = typename GET_PROP_TYPE(TypeTag, VtkOutputFields);
    VtkOutputFields::init(vtkWriter); //! Add model specific output fields
    vtkWriter.write(0.0);

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // run instationary non-linear simulation
    /////////////////////////////////////////////////////////////////////////////////////////////////

    //! set some check points for the time loop
    timeLoop->setPeriodicCheckPoint(tEnd/10.0);

    //! start the time loop
    timeLoop->start(); do
    {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

        Dune::Timer assembleTimer;
        assembler->assembleJacobianAndResidual(x);
        assembleTimer.stop();

        // solve the linear system A(xOld-xNew) = r
        Dune::Timer solveTimer;
        SolutionVector xDelta(x);
        linearSolver->solve(*A, xDelta, *r);
        solveTimer.stop();

        // update solution and grid variables
        Dune::Timer updateTimer;
        x -= xDelta;
        gridVariables->update(x);
        updateTimer.stop();

        // statistics
        const auto elapsedTot = assembleTimer.elapsed() + solveTimer.elapsed() + updateTimer.elapsed();
        std::cout << "Assemble/solve/update time: "
                  <<  assembleTimer.elapsed() << "(" << 100*assembleTimer.elapsed()/elapsedTot << "%)/"
                  <<  solveTimer.elapsed() << "(" << 100*solveTimer.elapsed()/elapsedTot << "%)/"
                  <<  updateTimer.elapsed() << "(" << 100*updateTimer.elapsed()/elapsedTot << "%)"
                  <<  std::endl;

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output on check points
        if (timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(dt);

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    //! print dumux end message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;

}
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
