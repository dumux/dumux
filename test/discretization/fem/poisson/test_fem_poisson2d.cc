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
 * \brief Test solving a 2d poisson problem with the finite element method.
 *        This helps identifying the minimum requirements for solving a problem
 *        using the finite element framework.
 */
#include <config.h>
#include <iostream>

#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/functions/gridfunctions/analyticgridviewfunction.hh>

#include <dumux/io/grid/gridmanager.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/discretization/fem/fegridgeometry.hh>
#include <dumux/discretization/fem/fegridvariables.hh>

#include "problem.hh"
#include "localresidual.hh"


int main (int argc, char *argv[]) try
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // initialize parameter tree
    Parameters::init(argc, argv);

    using Grid = Dune::YaspGrid<2>;
    using GridManager = GridManager<Grid>;
    GridManager gridManager;
    gridManager.init();

    using GridView = typename Grid::LeafGridView;
    const GridView& gridView = gridManager.grid().leafGridView();

    using Scalar = double;
    static constexpr int basisOrder = 1;

    // make the grid geometry
    using FEBasis = Dune::Functions::LagrangeBasis<GridView, basisOrder, Scalar>;
    using GridGeometry = FEGridGeometry<FEBasis>;
    auto feBasis = std::make_shared<FEBasis>(gridView);
    auto gridGeometry = std::make_shared<GridGeometry>(feBasis);
    gridGeometry->update();

    // make the problem
    using Tensor = TENSORTYPE; // Defined in CMakeLists.txt
    using Problem = FEPoissonProblem<Tensor, Scalar, GridGeometry>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // make the grid variables
    using PrimaryVariables = typename ProblemTraits<Problem>::PrimaryVariables;
    using SolutionVector = Dune::BlockVector<PrimaryVariables>;
    using GridVariables = FEGridVariables<Problem, SolutionVector>;
    auto gridVariables = std::make_shared<GridVariables>(problem);

    // initialize grid variables with initial solution
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;
    gridVariables->init(x);

    // make the assembler
    using LocalResidual = FEPoissonLocalResidual<GridVariables>;
    using Assembler = Assembler<GridVariables, LocalResidual, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(gridVariables);

    // non-linear solver (jacobian is inexact due to numeric differentiation)
    using LinearSolver = ILU0BiCGSTABBackend;
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    auto linearSolver = std::make_shared<LinearSolver>();
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver);

    // solve the system
    newtonSolver->solve(x);

    // Write out the numerical and the exact solutions
    auto evalExact = [problem] (const auto& pos) { return problem->exact(pos); };
    auto exactGF = Dune::Functions::makeAnalyticGridViewFunction(evalExact, gridGeometry->gridView());
    auto numericGF = Dune::Functions::makeDiscreteGlobalBasisFunction<PrimaryVariables>(gridGeometry->feBasis(), x);

    Dune::VTKWriter<GridView> vtkWriter(gridGeometry->gridView());
    vtkWriter.addVertexData(exactGF, Dune::VTK::FieldInfo({"x_numeric", Dune::VTK::FieldInfo::Type::scalar, 1}));
    vtkWriter.addVertexData(numericGF, Dune::VTK::FieldInfo({"x_exact", Dune::VTK::FieldInfo::Type::scalar, 1}));
    vtkWriter.write(problem->name());

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dune::Exception &e) {
    std::cout << e << std::endl;
    return 1;
}
