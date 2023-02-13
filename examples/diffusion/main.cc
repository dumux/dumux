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

#include <type_traits>
#include <random>

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/random.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/discretization/box.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/pdesolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include "model.hh"

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

namespace Dumux::Properties::TTag {

struct DiffusionTest
{
    using InheritsFrom = std::tuple<DiffusionModel, BoxModel>;

    using Scalar = double;
    using Grid = Dune::YaspGrid<2>;

    template<class TypeTag>
    using Problem = DiffusionTestProblem<TypeTag>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

} // end namespace Dumux::Properties

// initialize the solution with a random field
// and make sure this works in parallel as well
template<class SolutionVector, class GridGeometry>
SolutionVector createInitialSolution(const GridGeometry& gg)
{
    SolutionVector sol(gg.numDofs());

    // generate random number and add processor offset
    std::mt19937 gen(0); // seed is 0 for deterministic results
    Dumux::SimpleUniformDistribution<> dis(0.0, 1.0);
    // the rank is zero for sequential runs
    const auto rank = gg.gridView().comm().rank();
    for (int n = 0; n < sol.size(); ++n)
        sol[n] = dis(gen) + rank;

    // take the value of the processor with the minimum rank
    if (gg.gridView().comm().size() > 1)
    {
        Dumux::VectorCommDataHandleMin<
            typename GridGeometry::VertexMapper,
            SolutionVector,
            GridGeometry::GridView::dimension
        > minHandle(gg.vertexMapper(), sol);
        gg.gridView().communicate(minHandle, Dune::All_All_Interface, Dune::ForwardCommunication);

        // remove processor offset
        for (int n = 0; n < sol.size(); ++n)
            sol[n][0] -= std::floor(sol[n][0]);
    }

    return sol;
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // initialize parameter tree including command line arguments
    Parameters::init(argc, argv);

    // We specify an alias for the model type tag
    // We will configure the assembler with this type tag that
    // we specialized all these properties for above and in the model definition
    using TypeTag = Properties::TTag::DiffusionTest;

    // Moreover, these properties can be extracted to yield type information
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    // initialize the grid
    GridManager<Grid> gridManager;
    gridManager.init();

    // create the finite volume grid geometry from the (leaf) grid view
    // the problem for the boundary conditions, a solution vector and gridvariables
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);
    auto sol = createInitialSolution<SolutionVector>(*gridGeometry);
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // initialize the vtk output and write out the initial concentration field
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.priVar(0); }, "c");
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(
        getParam<Scalar>("TimeLoop.TStart", 0.0),
        getParam<Scalar>("TimeLoop.Dt"),
        getParam<Scalar>("TimeLoop.TEnd")
    );

    // Choose the type of assembler, linear solver and PDE solver
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = SSORCGIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::LinearPDESolver<Assembler, LinearSolver>;

    // Construct assembler, linear solver and PDE solver
    auto oldSol = sol; // copy the vector to store state of previous time step
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, oldSol);
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    Solver solver(assembler, linearSolver);

    // Assemble the system matrix once and then tell the solver to reuse in every time step
    // Since we have a linear problem, the matrix is not going to change
    assembler->setLinearSystem();
    assembler->assembleJacobian(sol);
    solver.reuseMatrix();

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

        // write VTK output (writes out the concentration field)
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

    } while (!timeLoop->finished());

    timeLoop->finalize(gridGeometry->gridView().comm());

    return 0;
}
