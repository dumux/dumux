// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

// In the file `main.cc`, we use the previously defined model to
// setup the simulation. The setup consists of four steps:
// 1. Define the problem setting boundary conditions
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
#include <dumux/io/chrono.hh>

// VTK output functionality
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

// the discretization method
#include <dumux/discretization/box.hh>

// assembly and solver
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>
#include <dumux/experimental/assembly/multistagefvassembler.hh>

// [[/codeblock]]
//
// Finally, we include the model defined in Part 1.
#include "model.hh"
// [[/details]]

//
// ## 1. The problem class
//
// The problem class implements the boundary conditions. It also provides
// an interface that is used by the local residual (see Part 1) to obtain the
// wave speed parameter. The value is read from the parameter configuration tree.
// [[content]]
namespace Dumux {
template<class TypeTag>
class WaveEquationProblem : public FVProblem<TypeTag>
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
    // In the constructor, we read the wave speed constant from the
    // parameter tree (which is initialized with the content of `params.input`).
public:
    WaveEquationProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        waveSpeed_ = Chrono::toSeconds(getParam("Problem.WaveSpeed")).count();
    }

    // As boundary conditions, we specify Neumann everywhere. This means, we
    // have to prescribe a flux at each boundary sub control volume face
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    // We prescribe zero flux over all of the boundary
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    // The wave speed interface is used in the local residual (see Part 1).
    // We can name this interface however we want as long as we adapt the calling site
    // in the `LocalResidual` class in `model.hh`.
    Scalar waveSpeed() const
    { return waveSpeed_; }

private:
    Scalar waveSpeed_;
};
} // end namespace Dumux
// [[/content]]

//
// ## 2. Property tag and specializations
//
// We create a new tag `WaveEquation` that inherits all properties
// specialized for the tag `WaveEquationModel` (we created this in Part 1)
// and the tag `BoxModel` which provides types relevant for the spatial
// discretization scheme (see [dumux/discretization/box.hh](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/discretization/box.hh)).
//
// Here we choose a short form of specializing properties. The property
// system also recognizes an alias (`using`) with the property name being
// a member of the specified type tag. Note that we could also use the same mechanism
// as in (Part 1), for example:
// ```code
// template<class TypeTag>
// struct Scalar<TypeTag, TTag::WaveEquation>
// { using type = double; };
//```
// which has the same effect as having an alias `Scalar = double;`
// as member of the type tag `WaveEquation`.
// This mechanism allows for a terser code expression.
// In case both definitions are present for the same type tag, the explicit
// specialization (long form) takes precedence.
//
// [[content]]
namespace Dumux::Properties::TTag {

struct WaveEquation
{
    using InheritsFrom = std::tuple<WaveEquationModel, BoxModel>;

    using Scalar = double;
    using Grid = Dune::YaspGrid<2>;

    template<class TypeTag>
    using Problem = WaveEquationProblem<TypeTag>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

} // end namespace Dumux::Properties::TTag
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
    // `./example_wave_equation -TimeLoop.TEnd 10s` will set the end time to $10$ seconds.
    // Command line arguments overwrite settings in the parameter file.
    Parameters::init(argc, argv);

    // We specify an alias for the model type tag.
    // We will configure the assembler with this type tag that
    // we specialized all these properties for above and in the model definition (Part 1).
    // We can extract type information through properties specialized for the type tag
    // using `GetPropType`.
    using TypeTag = Properties::TTag::WaveEquation;

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

    // initial solution (set all displacement dofs connected to the element
    // in the center of the domain to 1.0 to create a perturbation)
    SolutionVector sol(gridGeometry->numDofs());
    const auto gridCenter = 0.5*getParam<Dune::FieldVector<double, 2>>("Grid.UpperRight");
    const auto elementIndices = intersectingEntities(gridCenter, gridGeometry->boundingBoxTree());
    for (const auto eIdx : elementIndices)
    {
        const auto element = gridGeometry->element(eIdx);
        const auto fvGeometry = localView(*gridGeometry).bindElement(element);
        for (const auto& scv : scvs(fvGeometry))
            sol[scv.dofIndex()] = 1.0;
    }

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    // We initialize the VTK output module and write out the initial concentration field
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, sol, problem->name());
    vtkWriter.addVolumeVariable([](const auto& vv){ return vv.priVar(0); }, "c");
    vtkWriter.write(0.0);

    // We instantiate time loop using start and end time as well as
    // the time step size from the parameter tree (`params.input`)
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(
        Chrono::toSeconds(getParam("TimeLoop.TStart", "0")),
        Chrono::toSeconds(getParam("TimeLoop.Dt")),
        Chrono::toSeconds(getParam("TimeLoop.TEnd"))
    );

    auto timeSteppingScheme = []() -> std::shared_ptr<Experimental::MultiStageMethod<double>>
    {
        const auto timeStepScheme = getParam<std::string>("TimeLoop.Scheme", "DIRK3");
        if (timeStepScheme == "DIRK3")
            return std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<double>>();
        else if (timeStepScheme == "RK4")
            return std::make_shared<Experimental::MultiStage::RungeKuttaExplicitFourthOrder<double>>();
        else if (timeStepScheme == "ImplicitEuler")
            return std::make_shared<Experimental::MultiStage::ImplicitEuler<double>>();
        else if (timeStepScheme == "CrankNicolson")
            return std::make_shared<Experimental::MultiStage::Theta<double>>(0.5);
        else
            DUNE_THROW(Dumux::ParameterException,
                "Unknown time stepping scheme. Possible values are "
                << "DIRK3 or RK4 or ImplicitEuler or CrankNicolson");
    }();

    // the assembler with time loop for a transient problem
    auto oldSol = sol;
    using Assembler = Experimental::MultiStageFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeSteppingScheme, oldSol);

    // the linear solver
    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LAT>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    // the PDE solver
    using PDESolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    auto pdeSolver = std::make_shared<PDESolver>(assembler, linearSolver);

    using TimeStepper = Experimental::MultiStageTimeStepper<PDESolver>;
    TimeStepper timeStepper(pdeSolver, timeSteppingScheme);

    // The time loop is where most of the actual computations happen.
    // We assemble and solve the linear system of equations, update the solution,
    // write the solution to a VTK file and continue until the we reach the
    // final simulation time.
    //
    // [[codeblock]]
    timeLoop->start(); do
    {
       // linearize & solve
        timeStepper.step(sol, timeLoop->time(), timeLoop->timeStepSize());

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
