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
 * \ingroup OnePTests
 * \brief Convergence test with analytic solution
 */
#include <config.h>
#include <iostream>
#include <iomanip>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/assembly/fvassembler.hh>

#include "properties.hh"

template<class Problem, class SolutionVector>
void printL2Error(const Problem& problem, const SolutionVector& x)
{
    using namespace Dumux;
    using Scalar = double;

    Scalar l2error = 0.0;

    for (const auto& element : elements(problem.gridGeometry().gridView()))
    {
        auto fvGeometry = localView(problem.gridGeometry());
        fvGeometry.bindElement(element);

        for (auto&& scv : scvs(fvGeometry))
        {
            const auto dofIdx = scv.dofIndex();
            const Scalar delta = x[dofIdx] - problem.analyticalSolution(scv.dofPosition())[2/*pressureIdx*/];
            l2error += scv.volume()*(delta*delta);
        }
    }
    using std::sqrt;
    l2error = sqrt(l2error);

    const auto numDofs = problem.gridGeometry().numDofs();
    std::ostream tmp(std::cout.rdbuf());
    tmp << std::setprecision(8) << "** L2 error (abs) for "
            << std::setw(6) << numDofs << " dofs "
            << std::scientific
            << "L2 error = " << l2error
            << std::endl;

    // write the norm into a log file
    std::ofstream logFile;
    logFile.open(problem.name() + ".log", std::ios::app);
    logFile << "[ConvergenceTest] L2(p) = " << l2error << std::endl;
    logFile.close();
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TYPETAG;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // parse command line arguments
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using VelocityOutput = GetPropType<TypeTag, Properties::VelocityOutput>;
    vtkWriter.addVelocityOutput(std::make_shared<VelocityOutput>(*gridVariables));
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    // create assembler & linear solver
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // linearize & solve
    nonLinearSolver.solve(x);

    // output result to vtk
    vtkWriter.write(1.0);

    printL2Error(*problem, x);

    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;

}
