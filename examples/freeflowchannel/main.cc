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

 // ## The main file (`main.cc`)
 // [[content]]
 //
 // ### Included header files
 // [[details]] includes
 // [[exclude]]
 // Some generic includes.
#include <config.h>
#include <ctime>
#include <iostream>
// [[/exclude]]

// These are DUNE helper classes related to parallel computations and file I/O
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>

// The following headers include functionality related to property definition or retrieval, as well as
// the retrieval of input parameters specified in the input file or via the command line.
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

// The following files contain the non-linear Newton solver, the available linear solver backends and the assembler for the linear
// systems arising from the staggered-grid discretization.
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/assembly/diffmethod.hh> // analytic or numeric differentiation

// The following class provides a convenient way of writing of dumux simulation results to VTK format.
#include <dumux/io/staggeredvtkoutputmodule.hh>
// The gridmanager constructs a grid from the information in the input or grid file.
// Many different Dune grid implementations are supported, of which a list can be found
// in `gridmanager.hh`.
#include <dumux/io/grid/gridmanager.hh>

// This class contains functionality for additional flux output.
#include <dumux/freeflow/navierstokes/staggered/fluxoversurface.hh>


// In this header, a `TypeTag` is defined, which collects
// the properties that are required for the simulation.
// It also contains the actual problem with initial and boundary conditions.
// For detailed information, please have a look
// at the documentation provided therein.
#include "properties.hh"
// [[/details]]
//

// ### The main function
// We will now discuss the main program flow implemented within the `main` function.
// At the beginning of each program using Dune, an instance of `Dune::MPIHelper` has to
// be created. Moreover, we parse the run-time arguments from the command line and the
// input file:
// [[codeblock]]
int main(int argc, char** argv) try
{
    using namespace Dumux;

    // The Dune MPIHelper must be instantiated for each program using Dune
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);
    // [[/codeblock]]

    // We define a convenience alias for the type tag of the problem. The type
    // tag contains all the properties that are needed to define the model and the problem
    // setup. Throughout the main file, we will obtain types defined for this type tag
    // using the property system, i.e. with `GetPropType`.
    using TypeTag = Properties::TTag::ChannelExample;

    // #### Step 1: Create the grid
    // The `GridManager` class creates the grid from information given in the input file.
    // This can either be a grid file, or in the case of structured grids, one can specify the coordinates
    // of the corners of the grid and the number of cells to be used to discretize each spatial direction.
    // [[codeblock]]
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    // We compute on the leaf grid view.
    const auto& leafGridView = gridManager.grid().leafGridView();
    // [[/codeblock]]

    // #### Step 2: Setting up and solving the problem
    // First, a finite volume grid geometry is constructed from the grid that was created above.
    // This builds the sub-control volumes (scv) and sub-control volume faces (scvf) for each element
    // of the grid partition.
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // We now instantiate the problem, in which we define the boundary and initial conditions.
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // We set a solution vector which consist of two parts: one part (indexed by `cellCenterIdx`)
    // is for the pressure degrees of freedom (`dofs`) living in grid cell centers. Another part
    // (indexed by `faceIdx`) is for degrees of freedom defining the normal velocities on grid cell faces.
    // We initialize the solution vector by what was defined as the initial solution of the the problem.
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    x[GridGeometry::cellCenterIdx()].resize(gridGeometry->numCellCenterDofs());
    x[GridGeometry::faceIdx()].resize(gridGeometry->numFaceDofs());
    problem->applyInitialSolution(x);

    // The grid variables are used store variables (primary and secondary variables) on sub-control volumes and faces (volume and flux variables).
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // We then initialize the predefined model-specific output vtk output.
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    // [[details]] calculation of surface fluxes
    // We set up two surfaces over which fluxes are calculated.
    // We determine the extensions [xMin,xMax]x[yMin,yMax] of the physical domain.
    // The first surface (added by the first call of addSurface) shall be placed at the middle of the channel.
    // If we have an odd number of cells in x-direction, there would not be any cell faces
    // at the position of the surface (which is required for the flux calculation).
    // In this case, we add half a cell-width to the x-position in order to make sure that
    // the cell faces lie on the surface. This assumes a regular cartesian grid.
    // The second surface (second call of addSurface) is placed at the outlet of the channel.
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

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    const auto p0middle = GlobalPosition{planePosMiddleX + offsetX, yMin};
    const auto p1middle = GlobalPosition{planePosMiddleX + offsetX, yMax};
    flux.addSurface("middle", p0middle, p1middle);

    const auto p0outlet = GlobalPosition{xMax, yMin};
    const auto p1outlet = GlobalPosition{xMax, yMax};
    flux.addSurface("outlet", p0outlet, p1outlet);
    // [[/details]]

    // We create and initialize the assembler for the stationary problem.
    // This is where the Jacobian matrix for the Newton solver is assembled.
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    // We use UMFPack as direct linear solver within each Newton iteration.
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // This example considers a linear problem (incompressible Stokes flow), therefore
    // the non-linear Newton solver is not really necessary.
    // For sake of generality, we nevertheless use it here such that the example can be easily
    // changed to a non-linear problem by switching on the inertia terms in the input file or by choosing a compressible fluid.
    // In the following piece of code we instantiate the non-linear newton solver and let it solve
    // the problem.
    // [[codeblock]]
    // alias for and instantiation of the newton solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // Solve the (potentially non-linear) system.
    nonLinearSolver.solve(x);
    // [[/codeblock]]
    // In the following we calculate mass and volume fluxes over the planes specified above
    // (you have to click to unfold the code showing how to set up the surface fluxes).
    flux.calculateMassOrMoleFluxes();
    flux.calculateVolumeFluxes();

    // #### Final Output
    // We write the VTK output and print the mass/energy/volume fluxes over the planes.
    // We conclude by printing the dumux end message.
    // [[codeblock]]
    vtkWriter.write(1.0);

    if (GetPropType<TypeTag, Properties::ModelTraits>::enableEnergyBalance())
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
        Parameters::print();

    return 0;
} // end main
// [[/codeblock]]
// #### Exception handling
// In this part of the main file we catch and print possible exceptions that could
// occur during the simulation.
// [[details]] exception handler
// [[codeblock]]
// errors related to run-time parameters
catch (Dumux::ParameterException &e)
{
    std::cerr << std::endl << e << " ---> Abort!" << std::endl;
    return 1;
}
// errors related to the parsing of Dune grid files
catch (Dune::DGFException & e)
{
    std::cerr << "DGF exception thrown (" << e <<
                 "). Most likely, the DGF file name is wrong "
                 "or the DGF file is corrupted, "
                 "e.g. missing hash at end of file or wrong number (dimensions) of entries."
                 << " ---> Abort!" << std::endl;
    return 2;
}
// generic error handling with Dune::Exception
catch (Dune::Exception &e)
{
    std::cerr << "Dune reported error: " << e << " ---> Abort!" << std::endl;
    return 3;
}
// other exceptions
catch (...)
{
    std::cerr << "Unknown exception thrown! ---> Abort!" << std::endl;
    return 4;
}
// [[/codeblock]]
// [[/details]]
// [[/content]]
