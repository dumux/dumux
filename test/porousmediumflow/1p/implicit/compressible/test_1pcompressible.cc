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

#include <dumux/implicit/cellcentered/assembler.hh>

#include "problem.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = TTAG(IncompressibleTestProblem);

    // some aliases for better readability
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ParameterTree = typename GET_PROP(TypeTag, ParameterTree);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
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

    // the problem (boundary conditions)
    auto problem = std::make_shared<Problem>(leafGridView);

    // the solution vector
    auto x = std::make_shared<SolutionVector>(leafGridView.size(0));
    auto xold = x;

    // the grid variables
    auto gridVariables = std::make_shared<GridVariables>();
    gridVariables->init(*problem, *fvGridGeometry, *x);

    // // TEST
    // for (const auto& element : elements(leafGridView))
    // {
    //     auto fvGeometry = localView(*fvGridGeometry);
    //     fvGeometry.bindElement(element);

    //     auto elemVolVars = localView(gridVariables->curGridVolVars());
    //     elemVolVars.bindElement(element, fvGeometry, x);

    //     for (const auto& scv : scvs(fvGeometry))
    //     {
    //         const auto volVars = elemVolVars[scv];
    //         std::cout << volVars.pressure(0) << " ";
    //     }
    // }
    // std::cout << std::endl;

    // make assemble and attach linear system
    auto assembler = std::make_shared<CCImplicitAssembler<TypeTag>>(problem, fvGridGeometry, gridVariables);
    // auto A = std::make_shared<JacobianMatrix>();
    // auto r = std::make_shared<SolutionVector>();
    // assembler->setLinearSystem(A, r);

    // assemble the local jacobian and the residual
    // Dune::Timer timer; std::cout << "Assembling linear system ..." << std::flush;
    // assembler->assembleJacobianAndResidual(*x, *xold);
    // std::cout << " took " << timer.elapsed() << " seconds." << std::endl;

    // // print matrioc
    // // Dune::printmatrix(std::cout, *A, "", "");
    // // Dune::printvector(std::cout, *r, "", "");

    // we solve Ax = -r
    // (*r) *= -1.0;

    // // solve the linear system
    // timer.reset(); std::cout << "Solving linear system ..." << std::flush;
    auto linearSolver = std::make_shared<ILU0BiCGSTABBackend<TypeTag>>(*problem);
    // auto linearSolver = std::make_shared<UMFPackBackend<TypeTag>>(*problem);
    // linearSolver->solve(*A, *x, *r);
    // std::cout << " took " << timer.elapsed() << " seconds." << std::endl;

    NewtonMethod<TypeTag> nonLinearSolver;
    NewtonController newtonController(leafGridView.comm());
    nonLinearSolver.solve(newtonController, *assembler, *linearSolver, *x, *xold);

    // output result to vtk
    Dune::VTKWriter<GridView> vtkwriter(leafGridView);
    vtkwriter.addCellData(*x, "p");
    vtkwriter.write("test_1pincompressible");

    return 0;

} // end main
