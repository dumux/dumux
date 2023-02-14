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

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Problem //////// (often in a separate file problem.hh) ////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>

namespace Dumux {

template<class TypeTag>
class DiffusionTestProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
public:
    DiffusionTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        diffusionCoefficient_ = getParam<Scalar>("Problem.DiffusionCoefficient");
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return { 0.0 }; }

    Scalar diffusionCoefficient() const
    { return diffusionCoefficient_; }

private:
    Scalar diffusionCoefficient_;
};

} // end namespace Dumux

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Test case properties/traits /// (often in a separate file properties.hh) //
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <dune/grid/yaspgrid.hh>
#include <dumux/discretization/box.hh>
#include "model.hh"

namespace Dumux::Properties {

namespace TTag {
struct DiffusionTest { using InheritsFrom = std::tuple<DiffusionModel, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::DiffusionTest>
{ using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DiffusionTest>
{ using type = DiffusionTestProblem<TypeTag>; };

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DiffusionTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DiffusionTest>
{ static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DiffusionTest>
{ static constexpr bool value = true; };

} // end namespace Dumux::Properties

//////////////////////////////////////////////////////////////////////////////

#include <random>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/amgbackend.hh>
#include <dumux/linear/pdesolver.hh>
#include <dumux/assembly/fvassembler.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::DiffusionTest;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameter tree
    Parameters::init(argc, argv);

    // initialize the grid
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector sol(gridGeometry->numDofs());

    // create a random initial field
    std::mt19937 gen(0); // seed is 0 for deterministic results
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int n = 0; n < sol.size(); ++n)
        sol[n] = dis(gen);

    // take the minimum value at processor borders
    VectorCommDataHandleMin<typename GridGeometry::VertexMapper, SolutionVector, GridGeometry::GridView::dimension>
        minHandle(gridGeometry->vertexMapper(), sol);
    leafGridView.communicate(minHandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);

    // copy the vector to store state of previous time step
    auto oldSol = sol;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto dt = getParam<Scalar>("TimeLoop.Dt");

    // initialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.concentration(); }, "c");
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);

    // the assembler with time loop for a transient problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, oldSol);

    // the linear solver
    using LinearSolver = AMGBiCGSTABBackend<LinearSolverTraits<GridGeometry>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    // the solver
    using Solver = Dumux::LinearPDESolver<Assembler, LinearSolver>;
    Solver solver(assembler, linearSolver);

    // time loop
    timeLoop->start(); do
    {
        // assemble & solve
        solver.solve(sol);

        // make the new solution the old solution
        oldSol = sol;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write VTK output
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    return 0;
}
