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

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/partial.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "properties.hh"
#include "../../l2error.hh"

template<class Problem, class SolutionVector, class GridGeometry>
void printL2Error(const Problem& problem, const SolutionVector& x, const GridGeometry& gridGeometry)
{
    using namespace Dumux;
    using Scalar = double;
    using TypeTag = Properties::TTag::PipeFlow;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    using L2Error = NavierStokesTestL2Error<Scalar, ModelTraits, PrimaryVariables>;
    const auto [l2errorAbs, l2errorRel] = L2Error::calculateL2Error(*problem, x);
    const int numCellCenterDofs = gridGeometry->numCellCenterDofs();
    const int numFaceDofs = gridGeometry->numFaceDofs();
    std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                << std::setw(6) << numCellCenterDofs << " cc dofs and " << numFaceDofs << " face dofs (total: " << numCellCenterDofs + numFaceDofs << "): "
                << std::scientific
                << "L2(p) = " << l2errorAbs[Indices::pressureIdx] << " / " << l2errorRel[Indices::pressureIdx]
                << " , L2(vx) = " << l2errorAbs[Indices::velocityXIdx] << " / " << l2errorRel[Indices::velocityXIdx]
                << " , L2(vy) = " << l2errorAbs[Indices::velocityYIdx] << " / " << l2errorRel[Indices::velocityYIdx]
                << std::endl;

    // write the norm into a log file
    std::ofstream logFile;
    logFile.open(problem->name() + ".log", std::ios::app);
    logFile << "[ConvergenceTest] L2(p) = " << l2errorAbs[Indices::pressureIdx] << " L2(vx) = " << l2errorAbs[Indices::velocityXIdx] << " L2(vy) = " << l2errorAbs[Indices::velocityYIdx] << std::endl;
    logFile.close();
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // Define the sub problem type tags
    using TypeTag = Properties::TTag::PipeFlow;

    // try to create a grid (from the given grid file or the input file)
    Dumux::GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& gridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridView);
    gridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector sol;
    sol[GridGeometry::cellCenterIdx()].resize(gridGeometry->numCellCenterDofs());
    sol[GridGeometry::faceIdx()].resize(gridGeometry->numFaceDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    nonLinearSolver.solve(sol);
    vtkWriter.write(1.0);

    printL2Error(problem, sol, gridGeometry);

    // print dumux end message
    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;

} // end main
