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
 * \brief Test for the one-phase CC model
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/assembly/fvassembler.hh>

#include "problem.hh"

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

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry (and make it periodic)
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    bool periodic = gridGeometry->isPeriodic();
    periodic = gridGeometry->gridView().comm().max(periodic);
    if (!periodic)
        DUNE_THROW(Dune::GridError, "Your grid is not periodic. Maybe the grid manager doesn't support periodic boundaries.");

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);
    problem->computePointSourceMap();

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    // make assemble and attach linear system
    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    using JacobianMatrix = GetPropType<TypeTag, Properties::JacobianMatrix>;
    auto A = std::make_shared<JacobianMatrix>();
    auto r = std::make_shared<SolutionVector>();
    assembler->setLinearSystem(A, r);

    Dune::Timer timer;
    // assemble the local jacobian and the residual
    Dune::Timer assemblyTimer;
    if (mpiHelper.rank() == 0) std::cout << "Assembling linear system ..." << std::flush;
    assembler->assembleJacobianAndResidual(x);
    assemblyTimer.stop();
    if (mpiHelper.rank() == 0) std::cout << " took " << assemblyTimer.elapsed() << " seconds." << std::endl;

    // we solve Ax = -r to save update and copy
    (*r) *= -1.0;

    // solve the linear system
    Dune::Timer solverTimer;
    using LinearSolver = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    if (mpiHelper.rank() == 0) std::cout << "Solving linear system using " + linearSolver->name() + "..." << std::flush;
    linearSolver->solve(*A, x, *r);
    solverTimer.stop();
    if (mpiHelper.rank() == 0) std::cout << " took " << solverTimer.elapsed() << " seconds." << std::endl;

    // the grid variables need to be up to date for subsequent output
    Dune::Timer updateTimer; if (mpiHelper.rank() == 0) std::cout << "Updating variables ..." << std::flush;
    gridVariables->update(x);
    if (mpiHelper.rank() == 0) std::cout << " took " << updateTimer.elapsed() << std::endl;

    // output result to vtk
    vtkWriter.write(1.0);

    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    if (mpiHelper.rank() == 0)
        std::cout << "Simulation took " << timer.elapsed() << " seconds on "
                  << comm.size() << " processes.\n"
                  << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;

}
