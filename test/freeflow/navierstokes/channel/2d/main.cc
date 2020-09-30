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
 * \brief Channel flow test for the staggered grid (Navier-)Stokes model.
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
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/staggeredvtkoutputmodule.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include "problem.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::ChannelTest;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // check if we are about to restart a previously interrupted simulation
    Scalar restartTime = getParam<Scalar>("Restart.Time", 0);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    x[GridGeometry::cellCenterIdx()].resize(gridGeometry->numCellCenterDofs());
    x[GridGeometry::faceIdx()].resize(gridGeometry->numFaceDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(restartTime, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);
    problem->setTimeLoop(timeLoop);

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // initialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    StaggeredVtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields

    const bool isStationary = getParam<bool>("Problem.IsStationary", false);

    if (problem->hasAnalyticalSolution())
    {
        vtkWriter.addField(problem->getAnalyticalPressureSolution(), "pressureExact");
        vtkWriter.addField(problem->getAnalyticalVelocitySolution(), "velocityExact");
        vtkWriter.addFaceField(problem->getAnalyticalVelocitySolutionOnFace(), "faceVelocityExact");
    }

    vtkWriter.write(restartTime);

    // the assembler with time loop for instationary problem
    using Assembler = StaggeredFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = isStationary ? std::make_shared<Assembler>(problem, gridGeometry, gridVariables)
                                  : std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // set up two surfaces over which fluxes are calculated
    FluxOverSurface<GridVariables,
                    SolutionVector,
                    GetPropType<TypeTag, Properties::ModelTraits>,
                    GetPropType<TypeTag, Properties::LocalResidual>> flux(*gridVariables, x);
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    const Scalar xMin = gridGeometry->bBoxMin()[0];
    const Scalar xMax = gridGeometry->bBoxMax()[0];
    const Scalar yMin = gridGeometry->bBoxMin()[1];
    const Scalar yMax = gridGeometry->bBoxMax()[1];

    // The first surface shall be placed at the middle of the channel.
    // If we have an odd number of cells in x-direction, there would not be any cell faces
    // at the position of the surface (which is required for the flux calculation).
    // In this case, we add half a cell-width to the x-position in order to make sure that
    // the cell faces lie on the surface. This assumes a regular cartesian grid.
    const Scalar planePosMiddleX = xMin + 0.5*(xMax - xMin);
    int numCellsX = getParam<std::vector<int>>("Grid.Cells")[0];

    const unsigned int refinement = getParam<unsigned int>("Grid.Refinement", 0);

    numCellsX *= (1<<refinement);

    const Scalar offsetX = (numCellsX % 2 == 0) ? 0.0 : 0.5*((xMax - xMin) / numCellsX);

    const auto p0middle = GlobalPosition{planePosMiddleX + offsetX, yMin};
    const auto p1middle = GlobalPosition{planePosMiddleX + offsetX, yMax};
    flux.addSurface("middle", p0middle, p1middle);

    // The second surface is placed at the outlet of the channel.
    const auto p0outlet = GlobalPosition{xMax, yMin};
    const auto p1outlet = GlobalPosition{xMax, yMax};
    flux.addSurface("outlet", p0outlet, p1outlet);

    if (isStationary)
    {
        Dune::Timer timer;

        // solve the non-linear system with time step control
        nonLinearSolver.solve(x);

        if(problem->hasAnalyticalSolution())
            problem->printL2Error(x);

        // write vtk output
        vtkWriter.write(1.0);

        // calculate and print mass fluxes over the planes
        flux.calculateMassOrMoleFluxes();
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

        // calculate and print volume fluxes over the planes
        flux.calculateVolumeFluxes();
        std::cout << "volume flux at middle is: " << flux.netFlux("middle")[0] << std::endl;
        std::cout << "volume flux at outlet is: " << flux.netFlux("outlet")[0] << std::endl;

        timer.stop();
    }
    else
    {
        // time loop
        timeLoop->start(); do
        {
            // solve the non-linear system with time step control
            nonLinearSolver.solve(x, *timeLoop);

            // make the new solution the old solution
            xOld = x;
            gridVariables->advanceTimeStep();

            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();

            // write vtk output
            vtkWriter.write(timeLoop->time());

            // calculate and print mass fluxes over the planes
            flux.calculateMassOrMoleFluxes();
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

            // calculate and print volume fluxes over the planes
            flux.calculateVolumeFluxes();
            std::cout << "volume flux at middle is: " << flux.netFlux("middle")[0] << std::endl;
            std::cout << "volume flux at outlet is: " << flux.netFlux("outlet")[0] << std::endl;

            // report statistics of this time step
            timeLoop->reportTimeStep();

            // set new dt as suggested by newton solver
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        } while (!timeLoop->finished());

        timeLoop->finalize(leafGridView.comm());
    }

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
