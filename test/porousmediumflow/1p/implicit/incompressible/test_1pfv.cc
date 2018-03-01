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

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/parsolverbackend.hh>

#include <dumux/common/properties.hh>
#include <dumux/nonlinear/newtonmethod.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/valgrind.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/defaultusagemessage.hh>

#include <dumux/assembly/fvassembler.hh>

#include "problem.hh"

int main(int argc, char** argv) try
{
    using namespace Dumux;

    using TypeTag = TTAG(TYPETAG);

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

    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    GridCreator::makeGrid();
    GridCreator::loadBalance();

    // we compute on the leaf grid view
    const auto& leafGridView = GridCreator::grid().leafGridView();

    // create the finite volume grid geometry
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // the problem (boundary conditions)
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    // the solution vector
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    SolutionVector x(fvGridGeometry->numDofs());

    // the grid variables
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module
    VtkOutputModule<TypeTag> vtkWriter(*problem, *fvGridGeometry, *gridVariables, x, problem->name());
    using VtkOutputFields = typename GET_PROP_TYPE(TypeTag, VtkOutputFields);
    VtkOutputFields::init(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0.0);

    // make assemble and attach linear system
    using Assembler = FVAssembler<TypeTag, NUMDIFFMETHOD>;
    auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
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
    using LinearSolver = ParallelSSORCGBackend<TypeTag>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, fvGridGeometry->dofMapper());

    if (mpiHelper.rank() == 0) std::cout << "Solving linear system using " + linearSolver->name() + "..." << std::flush;
    bool converged = linearSolver->solve(*A, x, *r);
    bool convergedRemote = leafGridView.comm().max(converged);
    if (!convergedRemote)
        DUNE_THROW(NumericalProblem, "Linear solver did not converge!");
    solverTimer.stop();
    if (mpiHelper.rank() == 0) std::cout << " took " << solverTimer.elapsed() << " seconds." << std::endl;

    // the grid variables need to be up to date for subsequent output
    Dune::Timer updateTimer;
    if (mpiHelper.rank() == 0) std::cout << "Updating variables ..." << std::flush;
    gridVariables->update(x);
    if (mpiHelper.rank() == 0) std::cout << " took " << updateTimer.elapsed() << std::endl;

    // output result to vtk
    vtkWriter.write(1.0);

    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    if (mpiHelper.rank() == 0)
    {
        std::cout << "Simulation took " << timer.elapsed() << " seconds on "
                  << comm.size() << " processes.\n"
                  << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

        Parameters::print();
    }

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
