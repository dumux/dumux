// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief A test for the shallow water model (wet dam break).
 */
#include <config.h>


#include <array>
#include <ctime>
#include <iostream>
#include <cmath>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>

#include <test/freeflow/shallowwater/errors.hh>

#include "properties.hh"

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using TypeTag = Properties::TTag::TYPETAG;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

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

    // the problem (initial and boundary conditions)
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    // the solution vector
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // initialize the vtk output module
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;

    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables,x, problem->name());
    vtkWriter.addField(problem->getExactWaterDepth(), "exactWaterDepth");
    vtkWriter.addField(problem->getExactVelocityX(), "exactVelocityX");
    problem->updateAnalyticalSolution(x,*gridVariables,0.0);
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler with time loop for instationary problem
    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    // the linear solver
    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>,
                                               LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());

    // the domain is closed, so the water volume is conserved
    // compensated summation, so that the round-off of the sum over many cells does not mask the conservation
    auto waterVolume = [&](const SolutionVector& sol) {
        Scalar volume = 0.0, compensation = 0.0;
        for (const auto& element : elements(leafGridView, Dune::Partitions::interior))
        {
            const Scalar cellVolume = sol[gridGeometry->elementMapper().index(element)][0]*element.geometry().volume();
            const Scalar sum = volume + cellVolume;
            using std::abs;
            compensation += abs(volume) >= abs(cellVolume) ? (volume - sum) + cellVolume : (cellVolume - sum) + volume;
            volume = sum;
        }
        return leafGridView.comm().sum(volume + compensation);
    };
    const auto initialWaterVolume = waterVolume(x);

    // the non-linear solver
    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);

    //! set some check point at the end of the time loop
    timeLoop->setCheckPoint(tEnd);

    // time loop
    timeLoop->start(); do
    {
        nonLinearSolver.solve(x,*timeLoop);

        // update the analytical solution
        problem->updateAnalyticalSolution(x,*gridVariables,timeLoop->time()+timeLoop->timeStepSize());

        // make the new solution the old solution
        xOld = x;
        gridVariables->advanceTimeStep();

        // advance the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        if (timeLoop->isCheckPoint())
            vtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by newton controller
        timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));


    } while (!timeLoop->finished());

    timeLoop->finalize(leafGridView.comm());

    // compare with the exact solution of the Riemann problem at the final time
    problem->updateAnalyticalSolution(x, *gridVariables, timeLoop->time());
    const auto& exactWaterDepth = problem->getExactWaterDepth();
    const auto& exactVelocityX = problem->getExactVelocityX();
    const auto errors = ShallowWaterTest::computeErrors(*gridGeometry, x, [&](const auto& element){
        const auto eIdx = gridGeometry->elementMapper().index(element);
        return std::array<Scalar, 3>{{ exactWaterDepth[eIdx], exactVelocityX[eIdx], 0.0 }};
    });
    if (mpiHelper.rank() == 0)
        ShallowWaterTest::printErrors(errors, timeLoop->time());

    const auto relativeVolumeChange = std::abs(waterVolume(x) - initialWaterVolume)/initialWaterVolume;
    if (mpiHelper.rank() == 0)
        std::cout << "Relative change of the water volume: " << relativeVolumeChange << std::endl;
    if (relativeVolumeChange > getParam<Scalar>("Problem.VolumeTolerance", 1e-12))
        DUNE_THROW(Dune::Exception, "Water volume not conserved, relative change " << relativeVolumeChange);

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
