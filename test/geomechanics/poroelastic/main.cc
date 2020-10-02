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
 * \ingroup GeomechanicsTests
 * \brief Test for the poro-elastic model.
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/io/file/vtk.hh>

#include "problem.hh"

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/amgbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

// function to evaluate the element stresses
template< class StressType,
          class Problem,
          class GridVariables,
          class SolutionVector,
          class SigmaStorage >
void assembleElementStresses(SigmaStorage& sigmaStorage,
                             SigmaStorage& effSigmaStorage,
                             const Problem& problem,
                             const typename GridVariables::GridGeometry& gridGeometry,
                             const GridVariables& gridVariables,
                             const SolutionVector& x)
{
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        auto fvGeometry = localView(gridGeometry);
        auto elemVolVars = localView(gridVariables.curGridVolVars());

        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, x);

        // evaluate flux variables cache at cell center
        using FluxVarsCache = typename GridVariables::GridFluxVariablesCache::FluxVariablesCache;
        FluxVarsCache fluxVarsCache;
        fluxVarsCache.update(problem, element, fvGeometry, elemVolVars, element.geometry().center());

        // get lame parameters, the pressure and compute stress tensor
        const auto sigma = StressType::stressTensor(problem, element, fvGeometry, elemVolVars, fluxVarsCache);
        const auto effSigma = StressType::effectiveStressTensor(problem, element, fvGeometry, elemVolVars, fluxVarsCache);

        // pass values into storage container
        using GridGeometry = typename GridVariables::GridGeometry;
        for (int dir = 0; dir < GridGeometry::GridView::dimension; ++dir)
        {
            const auto eIdx = gridGeometry.elementMapper().index(element);
            sigmaStorage[dir][eIdx] = sigma[dir];
            effSigmaStorage[dir][eIdx] = effSigma[dir];
        }
    }
}

// main function
int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TestPoroElastic;

    // stop time for the entire computation
    Dune::Timer timer;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module and output fields
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    using VtkOutputModule = Dumux::VtkOutputModule<GridVariables, SolutionVector>;
    VtkOutputModule vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter);

    // also, add exact solution to the output
    SolutionVector xExact(gridGeometry->numDofs());
    for (const auto& v : vertices(leafGridView))
        xExact[ gridGeometry->vertexMapper().index(v) ] = problem->exactSolution(v.geometry().center());
    vtkWriter.addField(xExact, "u_exact");

    // Furthermore, write out element stress tensors
    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ForceVector = Dune::FieldVector< Scalar, dimWorld >;

    // containers to store sigma/effective sigma
    std::array< std::vector<ForceVector>, dim > sigmaStorage;
    std::array< std::vector<ForceVector>, dim > effSigmaStorage;

    const auto numCells = gridGeometry->gridView().size(0);
    std::for_each(sigmaStorage.begin(), sigmaStorage.end(), [numCells] (auto& sigma) { sigma.resize(numCells); });
    std::for_each(effSigmaStorage.begin(), effSigmaStorage.end(), [numCells] (auto& effSigma) { effSigma.resize(numCells); });

    for (int dir = 0; dir < dim; ++dir)
    {
        vtkWriter.addField(sigmaStorage[dir], "sigma_" + std::to_string(dir), VtkOutputModule::FieldType::element);
        vtkWriter.addField(effSigmaStorage[dir], "effSigma_" + std::to_string(dir), VtkOutputModule::FieldType::element);
    }

    // use convenience function to compute stresses
    using StressType = GetPropType<TypeTag, Properties::StressType>;
    assembleElementStresses<StressType>(sigmaStorage, effSigmaStorage, *problem, *gridGeometry, *gridVariables, x);

    // write initial solution
    vtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    // the linear solver
    using LinearSolver = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;
    auto linearSolver = std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // linearize & solve
    nonLinearSolver.solve(x);

    // the grid variables need to be up to date for subsequent output
    gridVariables->update(x);

    // write vtk output
    assembleElementStresses<StressType>(sigmaStorage, effSigmaStorage, *problem, *gridGeometry, *gridVariables, x);
    vtkWriter.write(1.0);

    // print time and say goodbye
    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    if (mpiHelper.rank() == 0)
        std::cout << "Simulation took " << timer.elapsed() << " seconds on "
                  << comm.size() << " processes.\n"
                  << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    // print parameters
    if (mpiHelper.rank() == 0)
        Parameters::print();

    return 0;
} // end main
