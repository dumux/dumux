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
 * \ingroup NavierStokesTests
 * \brief 3D Channel flow test for the staggered grid (Navier-)Stokes model
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/assembly/staggeredfvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/staggered/fluxoversurface.hh>
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "problem.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::ThreeDChannelTest;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // create a grid
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    Dumux::GridManager<Grid> gridManager;

#if HAVE_DUNE_SUBGRID && GRID_DIM == 3
    const bool isStaircaseGeometry = getParam<bool>("Problem.IsStaircaseGeometry", false);

    auto selector = [&](const auto& element)
    {
        if (!isStaircaseGeometry)
            return true;

        const Scalar deltaX = 0.003;
        const Scalar deltaZ = 0.000075;
        const Scalar deltaY = 0.0003;

        const Scalar eps = 1e-8;
        const auto globalPos = element.geometry().center();

        return globalPos[2] > (deltaZ/deltaX * globalPos[0] + deltaZ/deltaY * globalPos[1] - deltaZ + eps);
    };

    gridManager.init(selector, "Internal");
#else
    gridManager.init();
#endif

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    x[GridGeometry::cellCenterIdx()].resize(gridGeometry->numCellCenterDofs());
    x[GridGeometry::faceIdx()].resize(gridGeometry->numFaceDofs());

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // intialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // set up two planes over which fluxes are calculated
    FluxOverSurface<GridVariables,
                    SolutionVector,
                    GetPropType<TypeTag, Properties::ModelTraits>,
                    GetPropType<TypeTag, Properties::LocalResidual>> flux(*gridVariables, x);
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using GlobalPosition = Dune::FieldVector<Scalar, GridView::dimensionworld>;

    const Scalar xMin = gridGeometry->bBoxMin()[0];
    const Scalar xMax = gridGeometry->bBoxMax()[0];
    const Scalar yMin = gridGeometry->bBoxMin()[1];
    const Scalar yMax = gridGeometry->bBoxMax()[1];

#if GRID_DIM == 3
    const Scalar zMin = gridGeometry->bBoxMin()[2];
    const Scalar zMax = gridGeometry->bBoxMax()[2];
#endif

    // The first plane shall be placed at the middle of the channel.
    // If we have an odd number of cells in x-direction, there would not be any cell faces
    // at the postion of the plane (which is required for the flux calculation).
    // In this case, we add half a cell-width to the x-position in order to make sure that
    // the cell faces lie on the plane. This assumes a regular cartesian grid.
    // The second plane is placed at the outlet of the channel.
#if GRID_DIM == 3
    const auto p0inlet = GlobalPosition{xMin, yMin, zMin};
    const auto p1inlet = GlobalPosition{xMin, yMax, zMin};
    const auto p2inlet = GlobalPosition{xMin, yMin, zMax};
    const auto p3inlet = GlobalPosition{xMin, yMax, zMax};
    flux.addSurface("inlet", p0inlet, p1inlet, p2inlet, p3inlet);
#else
    const auto p0inlet = GlobalPosition{xMin, yMin};
    const auto p1inlet = GlobalPosition{xMin, yMax};
    flux.addSurface("inlet", p0inlet, p1inlet);
#endif

    const Scalar planePosMiddleX = xMin + 0.5*(xMax - xMin);
    int numCellsX = getParam<std::vector<int>>("Grid.Cells")[0];

    const unsigned int refinement = getParam<unsigned int>("Grid.Refinement", 0);

    numCellsX *= (1<<refinement);

    const Scalar offsetX = (numCellsX % 2 == 0) ? 0.0 : 0.5*((xMax - xMin) / numCellsX);

#if GRID_DIM == 3
    const auto p0middle = GlobalPosition{planePosMiddleX + offsetX, yMin, zMin};
    const auto p1middle = GlobalPosition{planePosMiddleX + offsetX, yMax, zMin};
    const auto p2middle = GlobalPosition{planePosMiddleX + offsetX, yMin, zMax};
    const auto p3middle = GlobalPosition{planePosMiddleX + offsetX, yMax, zMax};
    flux.addSurface("middle", p0middle, p1middle, p2middle, p3middle);
#else
    const auto p0middle = GlobalPosition{planePosMiddleX + offsetX, yMin};
    const auto p1middle = GlobalPosition{planePosMiddleX + offsetX, yMax};
flux.addSurface("middle", p0middle, p1middle);
#endif

    // The second plane is placed at the outlet of the channel.
#if GRID_DIM == 3
    const auto p0outlet = GlobalPosition{xMax, yMin, zMin};
    const auto p1outlet = GlobalPosition{xMax, yMax, zMin};
    const auto p2outlet = GlobalPosition{xMax, yMin, zMax};
    const auto p3outlet = GlobalPosition{xMax, yMax, zMax};
    flux.addSurface("outlet", p0outlet, p1outlet, p2outlet, p3outlet);
#else
    const auto p0outlet = GlobalPosition{xMax, yMin};
    const auto p1outlet = GlobalPosition{xMax, yMax};
    flux.addSurface("outlet", p0outlet, p1outlet);
#endif

    // linearize & solve
    Dune::Timer timer;
    nonLinearSolver.solve(x);

    // write vtk output
    vtkWriter.write(1.0);

    // calculate and print mass fluxes over the planes
    flux.calculateMassOrMoleFluxes();
    if(GetPropType<TypeTag, Properties::ModelTraits>::enableEnergyBalance())
    {
        std::cout << "mass / energy flux at inlet is: " << flux.netFlux("inlet") << std::endl;
        std::cout << "mass / energy flux at middle is: " << flux.netFlux("middle") << std::endl;
        std::cout << "mass / energy flux at outlet is: " << flux.netFlux("outlet") << std::endl;
    }
    else
    {
        std::cout << "mass flux at inlet is: " << flux.netFlux("inlet") << std::endl;
        std::cout << "mass flux at middle is: " << flux.netFlux("middle") << std::endl;
        std::cout << "mass flux at outlet is: " << flux.netFlux("outlet") << std::endl;
    }

    // calculate and print volume fluxes over the planes
    flux.calculateVolumeFluxes();
    std::cout << "volume flux at inlet is: " << flux.netFlux("inlet")[0] << std::endl;
    std::cout << "volume flux at middle is: " << flux.netFlux("middle")[0] << std::endl;
    std::cout << "volume flux at outlet is: " << flux.netFlux("outlet")[0] << std::endl;


    timer.stop();

    std::cout << "analyticalFlux: " << problem->analyticalFlux() << std::endl;

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";


    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
} // end main
