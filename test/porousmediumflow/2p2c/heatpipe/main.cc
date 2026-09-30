// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux-Lecture contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include "config.h"
#include <algorithm>
#include <ctime>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>

#include <dumux/discretization/method.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::HeatPipeTypeTag;

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////////////////////
    // try to create a grid (from the given grid file or the input file)
    /////////////////////////////////////////////////////////////////////

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto fvGridGeometry = std::make_shared<FVGridGeometry>(leafGridView);

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(fvGridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(fvGridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, fvGridGeometry);
    gridVariables->init(x);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");
    auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");

    // initialize the vtk output module
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter); //!< Add model specific output fields
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, fvGridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<FVGridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    // set some check points for the time loop
    timeLoop->setPeriodicCheckPoint(getParam<Scalar>("TimeLoop.EpisodeLength"));

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
        if (timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the newton solver
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    // compare the dry-out front, i.e. the first vertex where only the gas phase is present,
    // to the semi-analytical solution (see test_heatpipe_odesolver.cc)
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    Scalar dryOutPosition = fvGridGeometry->bBoxMax()[0];
    for (const auto& vertex : vertices(leafGridView))
    {
        const auto dofIdx = fvGridGeometry->vertexMapper().index(vertex);
        if (x[dofIdx].state() == Indices::secondPhaseOnly)
            dryOutPosition = std::min(dryOutPosition, vertex.geometry().center()[0]);
    }

    // numerical diffusion smears the front towards the heat source, so the simulated
    // front is expected to lie downstream of the semi-analytical one
    const auto referenceDryOutPosition = getParam<Scalar>("Problem.ReferenceDryOutPosition");
    const auto maxDeviation = getParam<Scalar>("Problem.MaxDryOutDeviation");
    const auto deviation = dryOutPosition - referenceDryOutPosition;
    std::cout << "Dry-out front at x = " << dryOutPosition << " m (semi-analytical: "
              << referenceDryOutPosition << " m, deviation: " << deviation << " m)" << std::endl;
    if (deviation < 0.0 || deviation > maxDeviation)
        DUNE_THROW(Dune::InvalidStateException, "Dry-out front deviates by " << deviation
                    << " m from the semi-analytical solution, which is outside the expected range [0, "
                    << maxDeviation << "] m");

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
