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

// We look now at the main file for the channel problem.
// ### Includes
// Necessary files are included.
//<details>
//  <summary>Toggle to expand details</summary>
//
// First, the Dune configuration file is include, the standard header file for C++ to get time and date information
// and another standard header for in- and output.
//
//<details>
//  <summary>Toggle to expand code (includes of problem file and of standard headers)</summary>
//
#include <config.h>

#include <ctime>
#include <iostream>
// </details>
//
// Dumux is based on DUNE, the Distributed and Unified Numerics Environment, which provides several grid managers and
// linear solvers. So we need some includes from that.
//<details>
//  <summary>Toggle to expand code (dune includes)</summary>
//
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>
// </details>
//
// In Dumux, a property system is used to specify the model. For this, different properties are defined containing
// type definitions, values and methods. All properties are declared in the file `properties.hh`.
// The file `parameters.hh` contains the parameter class, which manages the definition of input parameters by a default
// value, the inputfile or the command line.
// The file `dumuxmessage.hh` contains the class defining the start and end message of the simulation.
// The file `valgrind.hh` contains debugging funcionality.
//
// The file `seqsolverbackend.hh` contains the class, which defines the sequential linear solver backends.
// The nonlinear Newton's method is included, as well as the assembler, which assembles the linear systems for staggered-grid finite volume schemes.
// The containing class in the file `diffmethod.hh` defines the different differentiation methods used to compute the derivatives of the residual.
//
// We need the class in `staggeredvtkoutputmodule.hh` to simplify the writing of dumux simulation data to VTK format.
// The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the
// different supported grid managers.
//
// The following class contains functionality for additional flux output to the console.
//<details>
//  <summary>Toggle to expand code (dumux includes)</summary>
//
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/valgrind.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/staggered/fluxoversurface.hh>
#include <dumux/freeflow/navierstokes/model.hh>

#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/material/components/constant.hh>

#include "problem.hh"
// </details>
//
// </details>

//
// ### Setup basic properties for our simulation
// We setup the DuMux properties for our simulation (click [here](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course/blob/master/slides/dumux-course-properties.pdf) for DuMux course slides on the property system) within the namespace Properties, which is a sub-namespace of Dumux.
// 1. For every test problem, a new `TypeTag` has to be created, which is done within the namespace `TTag` (subnamespace of `Properties`). It inherits from the Navier-Stokes flow model and the staggered-grid discretization scheme.
// 2. The grid is chosen to be a two-dimensional YASP grid.
// 3. We set the `FluidSystem` to be a one-phase liquid with a single component. The class `Component::Constant` refers to a component with constant fluid properties (density, viscosity, ...) that can be set via the input file in the group `[0.Component]` where the number is the identifier given as template argument to the class template `Component::Constant`.
// 4. The problem class `ChannelExampleProblem`, which is forward declared before we enter `namespace Dumux` and defined later in this file, is defined to be the problem used in this test problem (charaterized by the TypeTag `ChannelExample`). The fluid system, which contains information about the properties such as density, viscosity or diffusion coefficient of the fluid we're simulating, is set to a constant one phase liquid.
// 5. We enable caching for the following classes (which stores values that were already calculated for later usage and thus results in higher memory usage but improved CPU speed): the grid volume variables, the grid flux variables, the finite volume grid geometry.
//
namespace Dumux::Properties {

namespace TTag {
struct ChannelExample { using InheritsFrom = std::tuple<NavierStokes, StaggeredFreeFlowModel>; };
} // namespace TTag

template<class TypeTag>
struct Grid<TypeTag, TTag::ChannelExample> { using type = Dune::YaspGrid<2>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::ChannelExample> { using type = Dumux::ChannelExampleProblem<TypeTag> ; };

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ChannelExample>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
};

template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ChannelExample> { static constexpr bool value = true; };
} // end namespace Dumux::Properties
//
//
// ### Beginning of the main function
// We begin the main function by making the type tag `ChannelExample`, that we defined in `problem.hh` for this test problem available here.
// Then we initializing the message passing interface (MPI), even if we do not plan to run the application in parallel. Finalizing of the MPI is done automatically on exit.
// We continue by printing the dumux start message and parsing the command line arguments and runtimeparameters from the input file in the init function.
//<details>
//  <summary>Toggle to expand code (beginning of main)</summary>
//
int main(int argc, char** argv) try
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::ChannelExample;

    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    Parameters::init(argc, argv);
    // </details>
    //
    // ### Set-up and solving of the problem
    //
    // A gridmanager tries to create the grid either from a grid file or the input file. You can learn more about grids in
    // Dumux in the [grid exercise](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux-course/-/tree/master/exercises/exercise-grids).
    // Then, we compute the leaf grid view, which is the overlap of all grid levels in hierarchical grids.
    // We create (`std::make_shared`) and initialize (`update`) the finite volume geometry to build up the subcontrolvolumes
    // and subcontrolvolume faces for each element of the grid partition.
    // In the problem, we define the boundary and initial conditions.
    // We set a solution vector which has a part (indexed by `cellCenterIdx`) for degrees of freedom (`dofs`) living
    // in grid cell centers - pressures - and a part (indexed by `faceIdx`) for degrees of freedom livin on grid cell faces.
    // We initialize the solution vector by what was defined as the initial solution in `problem.hh` (`applyInitialSolution`)
    // and then use the solution vector to intialize the `gridVariables`. Grid variables in general contain the s
    // primary variables (velocities, pressures) as well as secondary variables (density, viscosity, ...).
    // We then initialize the vtkoutput. Each model has a predefined model-specific output with relevant parameters
    // for that model. Here, it is pressure, velocity, density and process rank (relevant in the case of parallelisation).
    //<details>
    //  <summary>Toggle to expand code</summary>
    //
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    x[GridGeometry::cellCenterIdx()].resize(gridGeometry->numCellCenterDofs());
    x[GridGeometry::faceIdx()].resize(gridGeometry->numFaceDofs());
    problem->applyInitialSolution(x);

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);
    // </details>
    //
    // We set up two surfaces over which fluxes are calculated.
    // We determine the extensions [xMin,xMax]x[yMin,yMax] of the physical domain.
    // The first surface (added by the first call of addSurface) shall be placed at the middle of the channel.
    // If we have an odd number of cells in x-direction, there would not be any cell faces
    // at the position of the surface (which is required for the flux calculation).
    // In this case, we add half a cell-width to the x-position in order to make sure that
    // the cell faces lie on the surface. This assumes a regular cartesian grid.
    // The second surface (second call of addSurface) is placed at the outlet of the channel.
    //<details>
    //  <summary>Toggle to expand code (addition of surfaces)</summary>
    //
    FluxOverSurface<GridVariables,
                    SolutionVector,
                    GetPropType<TypeTag, Properties::ModelTraits>,
                    GetPropType<TypeTag, Properties::LocalResidual>> flux(*gridVariables, x);

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    const Scalar xMin = gridGeometry->bBoxMin()[0];
    const Scalar xMax = gridGeometry->bBoxMax()[0];
    const Scalar yMin = gridGeometry->bBoxMin()[1];
    const Scalar yMax = gridGeometry->bBoxMax()[1];

    const Scalar planePosMiddleX = xMin + 0.5*(xMax - xMin);
    int numCellsX = getParam<std::vector<int>>("Grid.Cells")[0];

    const unsigned int refinement = getParam<unsigned int>("Grid.Refinement", 0);
    numCellsX *= (1<<refinement);

    const Scalar offsetX = (numCellsX % 2 == 0) ? 0.0 : 0.5*((xMax - xMin) / numCellsX);

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    const auto p0middle = GlobalPosition{planePosMiddleX + offsetX, yMin};
    const auto p1middle = GlobalPosition{planePosMiddleX + offsetX, yMax};
    flux.addSurface("middle", p0middle, p1middle);

    const auto p0outlet = GlobalPosition{xMax, yMin};
    const auto p1outlet = GlobalPosition{xMax, yMax};
    flux.addSurface("outlet", p0outlet, p1outlet);
    // </details>
    //
    // The incompressible Stokes equation depends only linearly on the velocity, so we could use a linear solver to solve the problem.
    // Here, we use the show the more general case which would also work for incompressible fluids or the
    // Navier-Stokes equation. We use non-linear Newton solver for the solution.
    // In each step of the Newton solver, we assemble the respective linear system, including the jacobian matrix and the residual by the
    // `assembler`. The linear systems are here solved by the direct linear solver `UMFPack`.
    // As a postprocessing, we calculate mass and volume fluxes over the planes specified above.
    //<details>
    //  <summary>Toggle to expand code (assembly, solution process, postprocessing)</summary>
    //
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    nonLinearSolver.solve(x);

    flux.calculateMassOrMoleFluxes();
    flux.calculateVolumeFluxes();
    // </details>
    //
    // ### Final Output
    // We write the vtk output and print the mass/energy/volume fluxes over the planes.
    // We conclude by printing the dumux end message. After the end of the main function,
    // possibly catched error messages are printed.
    //<details>
    //  <summary>Toggle to expand code (final output)</summary>
    //
    vtkWriter.write(1.0);

    if(GetPropType<TypeTag, Properties::ModelTraits>::enableEnergyBalance())
    {
        std::cout << "mass / energy flux at middle is: " << flux.netFlux("middle") << std::endl;
        std::cout << "mass / energy flux at outlet is: " << flux.netFlux("outlet") << std::endl;
    }
    else
    {
        std::cout << "mass flux at middle is: " << flux.netFlux("middle") << std::endl;
        std::cout << "mass flux at outlet is: " << flux.netFlux("outlet") << std::endl;
    }

    std::cout << "volume flux at middle is: " << flux.netFlux("middle")[0] << std::endl;
    std::cout << "volume flux at outlet is: " << flux.netFlux("outlet")[0] << std::endl;

    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
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
// </details>
//
