// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

// # Cahn-Hilliard equation main program
//
// The file [`main.cc`](main.cc) contains the problem class `CahnHilliardTestProblem`,
// properties and traits specific to the test case as well as the `main` function.
// The setup consists of four steps:
// 1. Define the problem setting boundary conditions and the diffusion coefficient
// 2. Configure the property system reusing the model defined in Part 1
// 3. Define a function for setting the random initial condition
// 4. The main program defining all steps of the program
//
// __Table of contents__
//
// [TOC]
//
// We start in `main.cc` with the necessary header includes:
// [[details]] includes
#include <config.h>
#include <type_traits>
#include <random>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>

#include <dumux/discretization/box.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include <dune/grid/yaspgrid.hh>
#include <dumux/io/chrono.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include "model.hh"
// [[/details]]

//
// ## 1. The problem class
//
// The problem class implements boundary conditions and the source term.
// It also provides an interface that is used by the local residual (see Part 1) to obtain
// the mobility and surface tension. The values are read from the parameter configuration tree.
// [[content]]
namespace Dumux {
template<class TypeTag>
class CahnHilliardTestProblem : public FVProblem<TypeTag>
{
    // [[details]] boilerplate code
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    // [[/details]]
    // In the constructor, we read the parameter constants from the
    // parameter tree (which is initialized with the content of `params.input`).
public:
    CahnHilliardTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        mobility_ = getParam<Scalar>("Problem.Mobility");
        surfaceTension_ = getParam<Scalar>("Problem.SurfaceTension");
        energyScale_ = getParam<Scalar>("Problem.EnergyScale");
    }

    // In the `source` function, we implement the derivative of the free energy.
    // This demonstrates how parts of the local residual can be split into model specific
    // parts (see `CahnHilliardModelLocalResidual`) and parts that might change from scenario to scenario.
    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        NumEqVector values(0.0);
        const auto& c = elemVolVars[scv].concentration();
        values[Indices::chemicalPotentialEqIdx] = -energyScale_*2.0*c*(2.0*c*c - 3*c + 1);
        return values;
    }

    // We choose boundary flux (or Neumann) conditions for all equations on the entire boundary,
    // while specifying zero flux for both equations.
    // [[codeblock]]
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return { 0.0, 0.0 }; }
    // [[/codeblock]]

    // The parameters interfaces are used in the local residual (see Part 1).
    // We can name this interface however we want as long as we adapt the calling site
    // in the `CahnHilliardModelLocalResidual` class in `model.hh`.
    // [[codeblock]]
    Scalar mobility() const
    { return mobility_; }

    Scalar surfaceTension() const
    { return surfaceTension_; }

private:
    Scalar mobility_;
    Scalar surfaceTension_;
    Scalar energyScale_;
};
} // end namespace Dumux
// [[/codeblock]]
// [[/content]]

//
// ## 2. Property tag and specializations
//
// We create a new tag `DiffusionTest` that inherits all properties
// specialized for the tag `DiffusionModel` (we created this in Part 1)
// and the tag `BoxModel` which provides types relevant for the spatial
// discretization scheme (see [dumux/discretization/box.hh](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/discretization/box.hh)).
//
// [[content]]
namespace Dumux::Properties::TTag {
struct CahnHilliardTest
{
    using InheritsFrom = std::tuple<CahnHilliardModel, BoxModel>;

    using Grid = Dune::YaspGrid<2>;

    template<class TypeTag>
    using Problem = CahnHilliardTestProblem<TypeTag>;

    using EnableGridGeometryCache = std::true_type;
    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
};
} // end namespace Dumux::Properties
// [[/content]]
//
// ## 3. Creating the initial solution
//
// To initialize with a random field in parallel, where each processor
// creates its own random number sequence, we need to communicate the
// resulting values on the processor border and overlap.
// See [Diffusion equation example](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/diffusion/docs/main.md).
// for details. Here in addition, we need to provide a custom scatter operation
// that handles vector types. We only need to communicate the first entry (concentration).
//
// [[content]]
// [[codeblock]]
// functor for data communication with MPI
struct MinScatter
{
    template<class A, class B>
    static void apply(A& a, const B& b)
    { a[0] = std::min(a[0], b[0]); }
};

// create the random initial solution
template<class SolutionVector, class GridGeometry>
SolutionVector createInitialSolution(const GridGeometry& gg)
{
    SolutionVector sol(gg.numDofs());

    // Generate random number and add processor offset
    // For sequential runs `rank` always returns `0`.
    std::mt19937 gen(0); // seed is 0 for deterministic results
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (int n = 0; n < sol.size(); ++n)
    {
        sol[n][0] = 0.42 + 0.02*(0.5-dis(gen)) + gg.gridView().comm().rank();
        sol[n][1] = 0.0;
    }

    // We take the value of the processor with the minimum rank
    // and subtract the rank offset
    if (gg.gridView().comm().size() > 1)
    {
        Dumux::VectorCommDataHandle<
            typename GridGeometry::VertexMapper,
            SolutionVector,
            GridGeometry::GridView::dimension,
            MinScatter
        > minHandle(gg.vertexMapper(), sol);
        gg.gridView().communicate(minHandle, Dune::All_All_Interface, Dune::ForwardCommunication);

        // Remove processor offset
        for (int n = 0; n < sol.size(); ++n)
            sol[n][0] -= std::floor(sol[n][0]);
    }
    return sol;
}
// [[/codeblock]]
// [[/content]]
//
// ## 4. The main program
//
// The `main` function sets up the simulation framework, initializes runtime parameters,
// creates the grid and storage vectors
// for the variables, primary and secondary. It specifies and constructs and assembler, which
// assembles the discretized residual and system matrix (Jacobian of the model residual), as well as
// a linear solver. A Newton method is used to solve the nonlinear equations where in each Newton iteration
// the Jacobian and the residual needs to be reassembled and the resulting linear system is solved.
// The time loop controls the time stepping. Here, we let the Newton solver suggest an adaption of
// the time step size based on the convergence history of the nonlinear solver.
//
// [[content]]
int main(int argc, char** argv)
{
    using namespace Dumux;

    // We initialize MPI and/or multithreading backend. When not running
    // in parallel the MPI setup is skipped.
    Dumux::initialize(argc, argv);

    // Then we initialize parameter tree.
    Parameters::init(argc, argv);

    // We create an alias for the type tag for this problem
    // and extract type information through properties.
    using TypeTag = Properties::TTag::CahnHilliardTest;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    // We initialize the grid, grid geometry, problem, solution vector, and grid variables.
    GridManager<Grid> gridManager;
    gridManager.init();

    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);
    auto sol = createInitialSolution<SolutionVector>(*gridGeometry);
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // We initialize the VTK output module and write out the initial concentration and chemical potential
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.concentration(); }, "c");
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.chemicalPotential(); }, "mu");
    vtkWriter.write(0.0);

    // We instantiate time loop using start and end time as well as
    // the time step size from the parameter tree (`params.input`)
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(
        Chrono::toSeconds(getParam("TimeLoop.TStart", "0")),
        Chrono::toSeconds(getParam("TimeLoop.InitialTimeStepSize")),
        Chrono::toSeconds(getParam("TimeLoop.TEnd"))
    );

    // We set the maximum time step size allowed in the adaptive time stepping scheme.
    timeLoop->setMaxTimeStepSize(Chrono::toSeconds(getParam("TimeLoop.MaxTimeStepSize")));

    // Next, we choose the type of assembler, linear solver and PDE solver
    // and construct instances of these classes.
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = SSORBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::NewtonSolver<Assembler, LinearSolver>;

    auto oldSol = sol;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, oldSol);
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    Solver solver(assembler, linearSolver);

    // The time loop is where most of the actual computations happen.
    // We assemble and solve the linear system of equations, update the solution,
    // write the solution to a VTK file and continue until the we reach the
    // final simulation time.
    //
    // [[codeblock]]
    timeLoop->start(); do
    {
        // Assemble & solve
        // Passing the time loop enables the solver to repeat the time step
        // with a reduced time step size if the Newton solver does not converge.
        solver.solve(sol, *timeLoop);

        // make the new solution the old solution
        oldSol = sol;
        gridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write VTK output (concentration field and chemical potential)
        vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the Newton solver
        timeLoop->setTimeStepSize(solver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    // print the final report
    timeLoop->finalize(gridGeometry->gridView().comm());
    return 0;
}
// [[/codeblock]]
// [[/content]]
