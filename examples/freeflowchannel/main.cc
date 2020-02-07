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
//  <summary>Click to toggle details</summary>
//
// First, the configuration file is include, then the problem, followed by the standard header file for C++ to get time and date information
// and another standard header for in- and output.
//
//<details>
//  <summary>Click to toggle details</summary>
//
#include <config.h>

#include "problem.hh"

#include <ctime>
#include <iostream>
// </details>
//
// Dumux is based on DUNE, the Distributed and Unified Numerics Environment, which provides several grid managers and
// linear solvers. So we need some includes from that.
//<details>
//  <summary>Click to toggle details</summary>
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
// The file parameters.hh contains the parameter class, which manages the definition of input parameters by a default
// value, the inputfile or the command line.
// The file `dumuxmessage.hh` contains the class defining the start and end message of the simulation.
// The file valgrind.hh contains debugging funcionality.
//<details>
//  <summary>Click to toggle details</summary>
//
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/valgrind.hh>
// </details>
//
// The file seqsolverbackend.hh contains the class, which defines the sequential linear solver backends.
// The nonlinear Newton's method is included, as well as the assembler, which assembles the linear systems for staggered-grid finite volume schemes.
// The containing class in the file diffmethod.hh defines the different differentiation methods used to compute the derivatives of the residual.
//<details>
//  <summary>Click to toggle details</summary>
//
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
// </details>
//
// We need the class in staggeredvtkoutputmodule.hh to simplify the writing of dumux simulation data to VTK format.
// The gridmanager constructs a grid from the information in the input or grid file. There is a specification for the
// different supported grid managers.
//<details>
//  <summary>Click to toggle details</summary>
//
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager.hh>
// </details>
//
// The following class contains functionality for additional flux output to the console.
//<details>
//  <summary>Click to toggle details</summary>
//
#include <dumux/freeflow/navierstokes/staggered/fluxoversurface.hh>
// </details>
//
// </details>
//
// ### Beginning of the main function
// We begin the main function by defining the type tag for this problem, initializing MPI (finalizing is done automatically on exit),
// printing the dumux start message and parsing the command line arguments and input file in the init function.
//<details>
//  <summary>Click to toggle details</summary>
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
    // ### Create the grid
    // A gridmanager tries to create the grid either from a grid file or the input file.
    // Then, we compute on the leaf grid view.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    const auto& leafGridView = gridManager.grid().leafGridView();
    // </details>
    //
    // ### Set-up and solving of the problem
    //
    // We create and initialize the finite volume grid geometry, the problem, the linear system, including the jacobian matrix, the residual and the solution vector and the gridvariables.
    // #### Set-up
    //
    // We need the finite volume geometry to build up the subcontrolvolumes (scv) and subcontrolvolume faces (scvf) for each element of the grid partition.
    // In the problem, we define the boundary and initial conditions.
    // We initialize the solution vector
    // and then use the solution vector to intialize the `gridVariables`.
    // We then initialize the vtkoutput. Each model has a predefined model-specific output with relevant parameters
    // for that model.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
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
    //  <summary>Click to toggle details</summary>
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
    // </details>
    //
    // #### Assembling the linear system
    // We set the assembler.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    // </details>
    //
    // #### Solution
    // We set the linear and non-linear solver, solve the non-linear system and calculate mass and volume fluxes over the planes.
    //<details>
    //  <summary>Click to toggle details</summary>
    //
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
    // possible catched error messages are printed.
    //<details>
    //  <summary>Click to toggle details</summary>
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
