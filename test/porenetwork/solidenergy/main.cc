// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief heat conduction test for the pore network model (with solid properties)
 */
#include <config.h>

#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/initialize.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/io/grid/porenetwork/griddata.hh>

#include <dumux/porenetwork/common/pnmvtkoutputmodule.hh>
#include <dumux/porenetwork/common/boundaryflux.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    using SolidTypeTag = Properties::TTag::PNMSolidModel;

    // initialize MPI, finalize is done automatically on exit
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    ////////////////////////////////////////////////////////////
    // parse the command line arguments and input file
    ////////////////////////////////////////////////////////////

    // parse command line arguments
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    using GridManager = PoreNetwork::GridManager<3>;
    GridManager gridManager;
    gridManager.init();

    // we compute on the leaf grid view
    const auto solidLeafGridView = gridManager.grid().leafGridView();
    const auto solidGridData = gridManager.getGridData();

    // create the finite volume grid geometry
    using SolidGridGeometry = GetPropType<SolidTypeTag, Properties::GridGeometry>;
    auto solidGridGeometry = std::make_shared<SolidGridGeometry>(solidLeafGridView,*solidGridData);

    // the spatial parameters
    using SolidSpatialParams = GetPropType<SolidTypeTag, Properties::SpatialParams>;
    auto solidSpatialParams = std::make_shared<SolidSpatialParams>(solidGridGeometry); // from common/spatialparams constructor

    // the problem (boundary conditions)
    using SolidProblem = GetPropType<SolidTypeTag, Properties::Problem>;
    auto solidProblem = std::make_shared<SolidProblem>(solidGridGeometry, solidSpatialParams);

    // the solution vector
    using GridView = typename SolidGridGeometry::GridView;
    using SolutionVector = GetPropType<SolidTypeTag, Properties::SolutionVector>;
    SolutionVector sol(solidLeafGridView.size(GridView::dimension));;
    solidProblem->applyInitialSolution(sol);
    auto solOld = sol;

    // the grid variables
    using SolidGridVariables = GetPropType<SolidTypeTag, Properties::GridVariables>;
    auto solidGridVariables = std::make_shared<SolidGridVariables>(solidProblem, solidGridGeometry);
    solidGridVariables->init(sol);

    // get some time loop parameters
    using Scalar = GetPropType<SolidTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd", 1.0);
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize", 1.0);
    auto dt = getParam<Scalar>("TimeLoop.DtInitial", 1.0);

    // initialize the vtk output modules
    using SolidVtkWriter = PoreNetwork::VtkOutputModule<SolidGridVariables, GetPropType<SolidTypeTag, Properties::FluxVariables>, SolutionVector>;
    SolidVtkWriter solidVtkWriter(*solidGridVariables, sol, solidProblem->name());
    using IOFields = GetPropType<SolidTypeTag, Properties::IOFields>;
    IOFields::initOutputModule(solidVtkWriter);

    solidVtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    static const bool isStationary = getParam<bool>("Problem.IsStationary", true);

    // the assembler
    using Assembler = FVAssembler<SolidTypeTag, DiffMethod::numeric>;
    auto assembler = isStationary ? std::make_shared<Assembler>(solidProblem, solidGridGeometry, solidGridVariables) //stationary case
                                  : std::make_shared<Assembler>(solidProblem, solidGridGeometry, solidGridVariables, timeLoop, solOld); // transient case -> timeloop needed

    using LinearSolver =  UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    if(isStationary) //check if stationary or transient problem
    {
        // solve the non-linear system without time step control
        nonLinearSolver.solve(sol);

        // write vtk output
        solidVtkWriter.write(1);

    }else{ //solve transient problem
        // time loop
        timeLoop->start(); do
        {
            // solve the non-linear system with time step control
            nonLinearSolver.solve(sol, *timeLoop);

            // make the new solution the old solution
            solOld = sol;
            solidGridVariables->advanceTimeStep();

            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();

            // write vtk output
            solidVtkWriter.write(timeLoop->time());

            // report statistics of this time step
            timeLoop->reportTimeStep();

            // set new dt as suggested by newton solver
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        } while (!timeLoop->finished());

        timeLoop->finalize(solidLeafGridView.comm());
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
}
