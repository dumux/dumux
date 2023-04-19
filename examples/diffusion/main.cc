// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

// In the file `main.cc`, we use the previously defined model to
// setup the simulation. The setup consists of four steps:
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

// Include the header containing std::true_type and std::false_type
#include <type_traits>

// Include the headers useful for the random initial solution
#include <random>
#include <dumux/common/random.hh>

// As Dune grid interface implementation we will use
// the structured parallel grid manager YaspGrid
#include <dune/grid/yaspgrid.hh>

// Common includes for problem and main
// [[codeblock]]
#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>

// VTK output functionality
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

// the discretization method
#include <dumux/discretization/box.hh>

// assembly and solver
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/pdesolver.hh>
#include <dumux/assembly/fvassembler.hh>
// [[/codeblock]]
//
// Finally, we include the model defined in Part 1.
#include "model.hh"
// [[/details]]

//
// ## 1. The problem class
//
// The problem class implements the boundary conditions. It also provides
// an interface that is used by the local residual (see Part 1) to obtain the diffusion
// coefficient. The value is read from the parameter configuration tree.
// [[content]]
namespace Dumux {
template<class TypeTag>
class DiffusionTestProblem : public FVProblem<TypeTag>
{
    // [[details]] boilerplate code
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridGeometry::LocalView::Element::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    // [[/details]]
    // In the constructor, we read the diffusion coefficient constant from the
    // parameter tree (which is initialized with the content of `params.input`).
public:
    DiffusionTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        diffusionCoefficient_ = getParam<Scalar>("Problem.DiffusionCoefficient");
    }

    // As boundary conditions, we specify Neumann everywhere. This means, we
    // have to prescribe a flux at each boundary sub control volume face
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    // We prescribe zero flux over all of the boundary
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return { 0.0 }; }

    // The diffusion coefficient interface is used in the local residual (see Part 1).
    // We can name this interface however we want as long as we adapt the calling site
    // in the `LocalResidual` class in `model.hh`.
    Scalar diffusionCoefficient() const
    { return diffusionCoefficient_; }

private:
    Scalar diffusionCoefficient_;
};
} // end namespace Dumux
// [[/content]]

//
// ## 2. Property tag and specializations
//
// We create a new tag `DiffusionTest` that inherits all properties
// specialized for the tag `DiffusionModel` (we created this in Part 1)
// and the tag `BoxModel` which provides types relevant for the spatial
// discretization scheme (see [dumux/discretization/box.hh](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/discretization/box.hh)).
//
// Here we choose a short form of specializing properties. The property
// system also recognizes an alias (`using`) with the property name being
// a member of the specified type tag. Note that we could also use the same mechanism
// as in (Part 1), for example:
// ```code
// template<class TypeTag>
// struct Scalar<TypeTag, TTag::DiffusionTest>
// { using type = double; };
//```
// which has the same effect as having an alias `Scalar = double;`
// as member of the type tag `DiffusionTest`.
// This mechanism allows for a terser code expression.
// In case both definitions are present for the same type tag, the explicit
// specialization (long form) takes precedence.
//
// [[content]]
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

} // end namespace Dumux::Properties::TTag
// [[/content]]

//
// ## 3. Creating the initial solution
//
// We want to initialize the entries of the solution vector $c_{h,B}$
// with a random field such that $c_{h,B} \sim U(0,1)$. When running
// with MPI in parallel this requires communication. With just one
// processor (`gg.gridView().comm().size() == 1`), we can just iterate
// over all entries of the solution vector and assign a uniformly
// distributed random number. However, with multiple processes, the
// solution vector only contains a subset of the degrees of freedom
// living on the processor. Moreover, some degrees of freedom exist
// on multiple processes. Each processor generates their own random
// number which means we might generate different numbers for the
// same degree of freedom. Therefore, we communicate the number.
// The idea is that each process adds its rank number (processes are
// numbered starting from 0) to its value. Then we take the minimum
// over all processes which will take return the value of the processor
// with the lowest rank. Afterwards, we subtract the added rank offset.
//
// [[content]]
// [[codeblock]]
template<class SolutionVector, class GridGeometry>
SolutionVector createInitialSolution(const GridGeometry& gg)
{
    SolutionVector sol(gg.numDofs());

    // Generate random number and add processor offset
    // For sequential run, `rank` always returns `0`.
    std::mt19937 gen(0); // seed is 0 for deterministic results
    Dumux::SimpleUniformDistribution<> dis(0.0, 1.0);

    const auto rank = gg.gridView().comm().rank();
    for (int n = 0; n < sol.size(); ++n)
        sol[n] = dis(gen) + rank;

    // We take the value of the processor with the minimum rank
    // and subtract the rank offset
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
// [[/codeblock]]
// [[/content]]
//
// ## 4. The main program
//
// [[content]]
int main(int argc, char** argv)
{
    using namespace Dumux;

    // First, we initialize MPI and the multithreading backend.
    // This convenience function takes care that everything is setup in the right order and
    // has to be called for every Dumux simulation. `Dumux::initialize` also respects
    // the environment variable `DUMUX_NUM_THREADS` to restrict to amount of available cores
    // for multi-threaded code parts (for example the assembly).
    Dumux::initialize(argc, argv);

    // We initialize parameter tree including command line arguments.
    // This will, per default, read all parameters from the configuration file `params.input`
    // if such as file exists. Then it will look for command line arguments. For example
    // `./example_diffusion -TimeLoop.TEnd 10` will set the end time to $10$ seconds.
    // Command line arguments overwrite settings in the parameter file.
    Parameters::init(argc, argv);

    // We specify an alias for the model type tag.
    // We will configure the assembler with this type tag that
    // we specialized all these properties for above and in the model definition (Part 1).
    // We can extract type information through properties specialized for the type tag
    // using `GetPropType`.
    using TypeTag = Properties::TTag::DiffusionTest;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    // First, we initialize the grid. Grid parameters are taken from the input file
    // from the group `[Grid]` if it exists. You can also pass any other parameter
    // group (e.g. `"MyGroup"`) to `init` and then it will look in `[MyGroup.Grid]`.
    GridManager<Grid> gridManager;
    gridManager.init();

    // We, create the finite volume grid geometry from the (leaf) grid view,
    // the problem for the boundary conditions, a solution vector and a grid variables instance.
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);
    auto sol = createInitialSolution<SolutionVector>(*gridGeometry);
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // We initialize the VTK output module and write out the initial concentration field
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.priVar(0); }, "c");
    vtkWriter.write(0.0);

    // We instantiate time loop using start and end time as well as
    // the time step size from the parameter tree (`params.input`)
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(
        getParam<Scalar>("TimeLoop.TStart", 0.0),
        getParam<Scalar>("TimeLoop.Dt"),
        getParam<Scalar>("TimeLoop.TEnd")
    );

    // Next, we choose the type of assembler, linear solver and PDE solver
    // and construct instances of these classes.
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = SSORCGIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::LinearPDESolver<Assembler, LinearSolver>;

    auto oldSol = sol; // copy the vector to store state of previous time step
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, oldSol);
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    Solver solver(assembler, linearSolver);

    // We tell the assembler to create the system matrix and vector objects. Then we
    // assemble the system matrix once and then tell the solver to reuse in every time step.
    // Since we have a linear problem, the matrix is not going to change.
    assembler->setLinearSystem();
    assembler->assembleJacobian(sol);
    solver.reuseMatrix();

    // The time loop is where most of the actual computations happen.
    // We assemble and solve the linear system of equations, update the solution,
    // write the solution to a VTK file and continue until the we reach the
    // final simulation time.
    //
    // [[codeblock]]
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
// [[/codeblock]]
// [[/content]]
