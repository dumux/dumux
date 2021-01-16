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
 * \brief The solver of the single-phase convergence test
 */
#ifndef DUMUX_INCOMPRESSIBLE_ONEP_CONVERGENCETEST_SOLVER_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_CONVERGENCETEST_SOLVER_HH

#include <iostream>
#include <string>

#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/pdesolver.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include "problem.hh"

// return type of solveRefinementLevel()
// stores the grid geometry and the produced solution
template<class TypeTag>
struct SolutionStorage
{
    using Grid = Dumux::GetPropType<TypeTag, Dumux::Properties::Grid>;
    using GridGeometry = Dumux::GetPropType<TypeTag, Dumux::Properties::GridGeometry>;
    using SolutionVector = Dumux::GetPropType<TypeTag, Dumux::Properties::SolutionVector>;

public:
    Dumux::GridManager<Grid> gridManager;
    std::shared_ptr<GridGeometry> gridGeometry;
    std::shared_ptr<SolutionVector> solution;
};

/*!
 * \brief Solves the problem for a given number of cells per direction.
 * \param numCells the number of cells per direction to be used
 * \return returns an object of SolutionStorage, which carries
 *         the grid and the solution after solving the problem
 */
template<class TypeTag>
SolutionStorage<TypeTag> solveRefinementLevel(int numCells)
{
    using namespace Dumux;

    // adapt the parameter tree to carry given number of cells
    const auto c = std::to_string(numCells);
    Parameters::init([&c] (auto& tree) { tree["Grid.Cells"] = c + " " + c; });

    // the returned object of this function
    // we create the grid in there directly
    SolutionStorage<TypeTag> storage;
    auto& gridManager = storage.gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // start timer
    Dune::Timer timer;

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    auto x = std::make_shared<SolutionVector>(gridGeometry->numDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(*x);

    // create assembler & linear solver
    using Assembler = FVAssembler<TypeTag, DiffMethod::analytic>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // solver the linear problem
    LinearPDESolver<Assembler, LinearSolver> solver(assembler,  linearSolver);
    solver.solve(*x);

    // maybe output result to vtk
    if (getParam<bool>("IO.WriteVTK", false))
    {
        VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, *x, problem->name());
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        IOFields::initOutputModule(vtkWriter); // Add model specific output fields

        // add exact solution
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        std::vector<Scalar> exact(gridGeometry->numDofs());

        for (const auto& element : elements(gridGeometry->gridView()))
        {
            auto fvGeometry = localView(*gridGeometry);
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
                exact[scv.dofIndex()] = problem->exact(scv.dofPosition());
        }

        vtkWriter.addField(exact, "p_exact");
        vtkWriter.write(1.0);
    }

    // fill storage and return
    storage.gridGeometry = gridGeometry;
    storage.solution = x;
    return storage;
}

#endif
