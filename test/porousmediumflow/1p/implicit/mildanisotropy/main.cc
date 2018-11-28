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
 * \brief convergence test for the one-phase model
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/timer.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include "problem.hh"

//! computes the L2 norm of the error of a problem
template< class GridView, class Problem, class SolutionVector >
typename SolutionVector::field_type
computeL2Norm(const GridView& gridView,
              const Problem& problem,
              const SolutionVector& x)
{
    // container to store results in
    using Scalar = typename SolutionVector::field_type;

    // integration order
    int order = Dumux::getParam<int>("L2Norm.Order");

    // integrate error in each cell
    Scalar norm = 0.0;
    Scalar denominator = 0.0;
    for (const auto& element : elements(gridView))
    {
        // make discrete element solution
        const auto elemSol = elementSolution(element, x, problem.fvGridGeometry());

        // integrate the pressure error over the element
        const auto eg = element.geometry();
        const auto& rule = Dune::QuadratureRules<Scalar, GridView::dimension>::rule(eg.type(), order);
        for (auto&& qp : rule)
        {
            // integration point in global coordinates
            const auto ip = eg.global(qp.position());

            // exact and discrete solution at integration point
            const auto u = problem.exactSolution(ip);
            const auto uh = evalSolution(element, eg, elemSol, ip);

            norm += (uh - u)*(uh - u)*qp.weight()*eg.integrationElement(qp.position());
            denominator += u*u*qp.weight()*eg.integrationElement(qp.position());
        }
    }

    // take the root of the norm and scale it
    using std::sqrt;
    norm = std::sqrt(norm/denominator);
    return norm;
}

//! solve the problem
double solve(int argc, char** argv)
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
    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);
    fvGridGeometry->update();

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(fvGridGeometry->numDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x);

    // Instantiate the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());

    // container for the output of the exact solution
    std::vector<GetPropType<TypeTag, Properties::Scalar>> exactSolution;

    // write output if desired
    const bool writeVtk = getParam<bool>("IO.WriteVtk");
    if (writeVtk)
    {
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        IOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
        vtkWriter.write(0.0);

        exactSolution.resize(fvGridGeometry->numDofs());
        for (const auto& element : elements(fvGridGeometry->gridView()))
        {
            if (FVGridGeometry::discMethod == DiscretizationMethod::box)
                for (int i = 0; i < element.geometry().corners(); ++i)
                    exactSolution[ fvGridGeometry->vertexMapper().subIndex(element, i, Grid::dimension) ]
                                 = problem->exactSolution( element.template subEntity<Grid::dimension>(i).geometry().center() );
            else
                exactSolution[ fvGridGeometry->elementMapper().index(element) ] = problem->exactSolution( element.geometry().center() );
        }
        vtkWriter.addField(exactSolution, "exact solution");
    }

    // make assemble and attach linear system
    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables);
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
    using LinearSolver = SSORCGBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    if (mpiHelper.rank() == 0) std::cout << "Solving linear system using " + linearSolver->name() + "..." << std::flush;
    linearSolver->solve(*A, x, *r);
    solverTimer.stop();
    if (mpiHelper.rank() == 0) std::cout << " took " << solverTimer.elapsed() << " seconds." << std::endl;

    // the grid variables need to be up to date for subsequent output
    Dune::Timer updateTimer; std::cout << "Updating variables ..." << std::flush;
    gridVariables->update(x);
    updateTimer.elapsed(); std::cout << " took " << updateTimer.elapsed() << std::endl;

    // output result to vtk
    if (writeVtk) vtkWriter.write(1.0);

    // print final message
    timer.stop();
    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    if (mpiHelper.rank() == 0)
        std::cout << "Simulation took " << timer.elapsed() << " seconds on "
                  << comm.size() << " processes.\n"
                  << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    if (mpiHelper.rank() == 0)
        Parameters::print();

    // compute and return l2 error
    return computeL2Norm(leafGridView, *problem, x);
}

// the main function
int main(int argc, char** argv) try
{
    int numRefinements = 10;

    // container to store the l2 errors & deltaX
    std::vector<double> deltaX(numRefinements);
    std::vector<double> errors(numRefinements);

    for (unsigned int i = 0; i < numRefinements; ++i)
    {
        std::cout << "Solving refinement " << i+1 << std::endl;

        // Cells parameter for this refinement level
        const int numCells = 10*(i+1);
        std::string argKey = "-Grid.Cells";
        std::string argValue = std::to_string(numCells) + " " + std::to_string(numCells);

        // create new argv array
        char** rArray;
        rArray = new char*[argc+3];
        for (int i = 0; i < argc; i++)
        {
            rArray[i] = new char[ strlen(argv[i]) + 1 ];
            strcpy(rArray[i], argv[i]);
        }

        // copy the key and value
        rArray[argc] = new char[argKey.size()];
        rArray[argc+1] = new char[argValue.size()];
        strcpy(rArray[argc], argKey.c_str());
        strcpy(rArray[argc+1], argValue.c_str());

        // ensure null pointer after the args
        rArray[argc+2] = NULL;

        // solve with the additional arguments
        deltaX[i] = 1.0/numCells;
        errors[i] = solve(argc+2, rArray);

        // free the reserved memory
        for (int i = 0; i < argc+2; i++)
            delete [] rArray[i];
        delete [] rArray;
    }

    std::cout << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << "Convergence test finished, computed the following L2 errors:" << std::endl;
    for (int i = 0; i < numRefinements; ++i)
        std::cout << "dx = " << deltaX[i] << "\t -> \t" << "l2 = " << errors[i] << std::endl;

    using std::log;
    std::cout << std::endl;
    std::cout << "Convergence rate: " << (log(errors[numRefinements-1]) - log(errors[numRefinements-2]))/
                                         (log(deltaX[numRefinements-1]) - log(deltaX[numRefinements-2])) << std::endl;
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
