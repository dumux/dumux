// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief A test for the shallow water model (bowl).
 *
 * The time integration scheme is `TimeLoop.Scheme`: `ImplicitEuler` (default) or `BDF2` (second-order backward
 * differentiation formula with variable step sizes, see Experimental::MultiStep::BDF).
 */
#include <config.h>

#include <cmath>
#include <ctime>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/experimental/assembly/multistagefvassembler.hh>
#include <dumux/experimental/timestepping/multistepmethods.hh>
#include <dumux/experimental/timestepping/multisteptimestepper.hh>

#include <test/freeflow/shallowwater/errors.hh>

#include "properties.hh"

template<class SolutionVector, int i>
class SolutionComponent
{
    const SolutionVector& sol_;
public:
    SolutionComponent(const SolutionVector& sol) : sol_(sol) {}
    auto operator[] (std::size_t idx) const { return sol_[idx][i]; }
    auto size() const { return sol_.size(); }
};

//! Compute the analytical solution at time t
template<class SolutionVector, class Problem>
SolutionVector computeAnalyticalSolution(const double t, const Problem& problem)
{
    const auto& gg = problem.gridGeometry();
    SolutionVector exactSolution(gg.numDofs());
    for (const auto& element : elements(gg.gridView()))
    {
        const auto eIdx = gg.elementMapper().index(element);
        const auto& globalPos = element.geometry().center();
        exactSolution[eIdx] = problem.analyticalSolution(t, globalPos);
    }
    return exactSolution;
}

//! The linear solver for the given assembler
template<class Assembler, class GridView, class GridGeometry>
auto makeLinearSolver(const GridView& gridView, const GridGeometry& gridGeometry)
{
    using namespace Dumux;
    using LinearSolver = LINEARSOLVER<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    if constexpr (std::is_constructible_v<LinearSolver, const GridView&, decltype(gridGeometry.dofMapper())>)
        return std::make_shared<LinearSolver>(gridView, gridGeometry.dofMapper());
    else
        return std::make_shared<LinearSolver>();
}

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
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    // the grid variables
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    // get some time loop parameters
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    const auto tEnd = 2.0*problem->oscillationPeriodInSeconds();
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // initialize the vtk output module and the analytical solution
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables,x, problem->name());
    auto exactSolution = computeAnalyticalSolution<SolutionVector>(0.0, *problem);
    const auto exactWaterDepth = SolutionComponent<SolutionVector, 0>(exactSolution);
    const auto exactVelocityX = SolutionComponent<SolutionVector, 1>(exactSolution);
    const auto exactVelocityY = SolutionComponent<SolutionVector, 2>(exactSolution);
    vtkWriter.addField(exactWaterDepth, "exactWaterDepth");
    vtkWriter.addField(exactVelocityX, "exactVelocityX");
    vtkWriter.addField(exactVelocityY, "exactVelocityY");
    using IOFields = GetPropType<TypeTag, Properties::IOFields>;
    IOFields::initOutputModule(vtkWriter);
    vtkWriter.write(0.0);

    // instantiate time loop
    auto timeLoop = std::make_shared<CheckPointTimeLoop<Scalar>>(0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    const auto exactSolutionOfElement = [&](const auto& element){
        return exactSolution[gridGeometry->elementMapper().index(element)];
    };
    const auto printErrors = [&](const Scalar t){
        const auto errors = ShallowWaterTest::computeErrors(*gridGeometry, x, exactSolutionOfElement);
        if (mpiHelper.rank() == 0)
            ShallowWaterTest::printErrors(errors, t);
        return errors;
    };

    // the bowl is closed, so the water volume is conserved
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

    // the initial solution is exact except for the thin water film on the dry bed
    printErrors(0.0);

    // time loop with a function that solves one time step and the non-linear solver suggesting the step size
    const auto runTimeLoop = [&](auto&& solveTimeStep, auto& nonLinearSolver)
    {
        timeLoop->start(); do
        {
            solveTimeStep();

            // make the new solution the old solution
            xOld = x;
            gridVariables->advanceTimeStep();

            // advance the time loop to the next step
            timeLoop->advanceTimeStep();

            // update the analytical solution and print the current errors
            exactSolution = computeAnalyticalSolution<SolutionVector>(timeLoop->time(), *problem);
            printErrors(timeLoop->time());

            // write vtk output every 10 time steps
            if (!(timeLoop->timeStepIndex() % 10))
                vtkWriter.write(timeLoop->time());

            // report statistics of this time step
            timeLoop->reportTimeStep();

            // set new dt as suggested by newton controller
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));

        } while (!timeLoop->finished());
    };

    const auto scheme = getParam<std::string>("TimeLoop.Scheme", "ImplicitEuler");
    if (scheme == "ImplicitEuler")
    {
        using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);
        auto linearSolver = makeLinearSolver<Assembler>(leafGridView, *gridGeometry);
        using NewtonSolver = Dumux::NewtonSolver<Assembler, typename decltype(linearSolver)::element_type>;
        NewtonSolver nonLinearSolver(assembler, linearSolver);
        runTimeLoop([&]{ nonLinearSolver.solve(x, *timeLoop); }, nonLinearSolver);
    }
    else if (scheme == "BDF2")
    {
        auto method = std::make_shared<Experimental::MultiStep::BDF<Scalar>>(2);
        using Assembler = Experimental::MultiStageFVAssembler<TypeTag, DiffMethod::numeric>;
        auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, method, xOld);
        auto linearSolver = makeLinearSolver<Assembler>(leafGridView, *gridGeometry);
        using NewtonSolver = Dumux::NewtonSolver<Assembler, typename decltype(linearSolver)::element_type>;
        auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver);
        Experimental::MultiStepTimeStepper<NewtonSolver> timeStepper(nonLinearSolver, method);

        // the first step and steps after a step size increase by more than a factor of two are taken with BDF1
        std::size_t numFirstOrderSteps = 0;
        runTimeLoop([&]{
            timeStepper.step(x, *timeLoop);
            numFirstOrderSteps += timeStepper.lastNumSteps() < 2;
        }, *nonLinearSolver);
        if (mpiHelper.rank() == 0)
            std::cout << "Time steps taken with BDF1: " << numFirstOrderSteps << std::endl;
    }
    else
        DUNE_THROW(ParameterException, "Unknown time integration scheme " << scheme << " (use ImplicitEuler or BDF2)");

    timeLoop->finalize(leafGridView.comm());

    const auto relativeVolumeChange = std::abs(waterVolume(x) - initialWaterVolume)/initialWaterVolume;
    if (mpiHelper.rank() == 0)
        std::cout << "Relative change of the water volume: " << relativeVolumeChange << std::endl;
    if (relativeVolumeChange > getParam<Scalar>("Problem.VolumeTolerance", 1e-12))
        DUNE_THROW(Dune::Exception, "Water volume not conserved, relative change " << relativeVolumeChange);

    const auto finalErrors = ShallowWaterTest::computeErrors(*gridGeometry, x, exactSolutionOfElement);
    const auto errorTolerance = getParam<Scalar>("Problem.ErrorTolerance", -1.0);
    if (errorTolerance > 0.0 && !(finalErrors.l1[0] < errorTolerance))
        DUNE_THROW(Dune::Exception, "L1 error of the water depth " << finalErrors.l1[0] << " exceeds " << errorTolerance);

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
