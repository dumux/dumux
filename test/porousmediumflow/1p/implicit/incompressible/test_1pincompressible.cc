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

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/common/propertysystem.hh>
#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>
#include <dumux/common/parameterparser.hh>

#include <dumux/assembly/ccassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include "problem.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = TTAG(IncompressibleTestProblem);

    // some aliases for better readability
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ParameterTree = typename GET_PROP(TypeTag, ParameterTree);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

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

    // the problem (boundary conditions)
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    // the solution vector
    auto x = std::make_shared<SolutionVector>(leafGridView.size(0));

    // the grid variables
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(*x);

    // make assemble and attach linear system
    auto assembler = std::make_shared<CCAssembler<TypeTag, DiffMethod::numeric>>(problem, fvGridGeometry, gridVariables);
    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<SolutionVector>();
    assembler->setLinearSystem(A, r);

    Dune::Timer timer;
    // assemble the local jacobian and the residual
    Dune::Timer assemblyTimer; std::cout << "Assembling linear system ..." << std::flush;
    assembler->assembleJacobianAndResidual(*x);
    assemblyTimer.stop(); std::cout << " took " << assemblyTimer.elapsed() << " seconds." << std::endl;

    // we solve Ax = -r
    (*r) *= -1.0;

    // // solve the linear system
    Dune::Timer solverTimer; std::cout << "Solving linear system ..." << std::flush;
    auto linearSolver = std::make_shared<ILU0BiCGSTABBackend<TypeTag>>(*problem);
    linearSolver->solve(*A, *x, *r);
    solverTimer.stop(); std::cout << " took " << solverTimer.elapsed() << " seconds." << std::endl;

    // output result to vtk
    Dune::VTKWriter<GridView> vtkwriter(leafGridView);
    vtkwriter.addCellData(*x, "p");
    vtkwriter.write("test_1pincompressible");

    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    return 0;

} // end main
