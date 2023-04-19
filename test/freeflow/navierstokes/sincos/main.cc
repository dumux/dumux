// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the instationary staggered grid Navier-Stokes model
 *        with analytical solution.
 */

#include <config.h>

#include <ctime>
#include <iostream>
#include <type_traits>
#include <tuple>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/vtkoutputmodule.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/istlsolverfactorybackend.hh>

#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/freeflow/navierstokes/velocityoutput.hh>
#include <test/freeflow/navierstokes/analyticalsolutionvectors.hh>
#include <test/freeflow/navierstokes/errors.hh>

#include "properties.hh"

template<class MomentumProblem>
auto createSource(const MomentumProblem& momentumProblem)
{
    using Scalar = double;
    using Indices = typename MomentumProblem::Indices;

    const auto& gridGeometry = momentumProblem.gridGeometry();
    std::array<std::vector<Scalar>, 2> source;

    for (auto& component : source)
        component.resize(gridGeometry.gridView().size(0));

    for (const auto& element : elements(gridGeometry.gridView()))
    {
        const auto& center = element.geometry().center();
        const auto eIdx = gridGeometry.elementMapper().index(element);
        auto sourceAtPosVal = momentumProblem.sourceAtPos(center);
        source[Indices::momentumXBalanceIdx][eIdx] = sourceAtPosVal[Indices::momentumXBalanceIdx];
        source[Indices::momentumYBalanceIdx][eIdx] = sourceAtPosVal[Indices::momentumYBalanceIdx];
    }

    return source;
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::SincosTestMomentum;
    using MassTypeTag = Properties::TTag::SincosTestMass;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<GetPropType<MomentumTypeTag, Properties::Grid>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

    // the coupling manager
    using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>();

    // get some time loop parameters
    using Scalar = GetPropType<MomentumTypeTag, Properties::Scalar>;
    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto maxDt = getParam<Scalar>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<Scalar>("TimeLoop.DtInitial");

    // instantiate time loop
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the problem (initial and boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);
    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
    constexpr auto massIdx = CouplingManager::freeFlowMassIndex;
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    momentumProblem->applyInitialSolution(x[momentumIdx]);
    massProblem->applyInitialSolution(x[massIdx]);
    auto xOld = x;

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);
    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    const bool isStationary = getParam<bool>("Problem.IsStationary");
    if (isStationary)
        couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    else
        couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x, xOld);

    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    // initialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    NavierStokesTest::AnalyticalSolutionVectors analyticalSolVectors(momentumProblem, massProblem);
    vtkWriter.addField(analyticalSolVectors.analyticalPressureSolution(), "pressureExact");
    vtkWriter.addField(analyticalSolVectors.analyticalVelocitySolution(), "velocityExact");
    const auto source = createSource(*momentumProblem);
    const auto& sourceX = source[MomentumProblem::Indices::momentumXBalanceIdx];
    const auto& sourceY = source[MomentumProblem::Indices::momentumYBalanceIdx];
    vtkWriter.addField(sourceX, "sourceX");
    vtkWriter.addField(sourceY, "sourceY");
    vtkWriter.write(0.0);

    // the assembler with time loop for instationary problem
    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = isStationary ?
        std::make_shared<Assembler>(
            std::make_tuple(momentumProblem, massProblem),
            std::make_tuple(momentumGridGeometry, massGridGeometry),
            std::make_tuple(momentumGridVariables, massGridVariables),
            couplingManager
        )
        :
        std::make_shared<Assembler>(
            std::make_tuple(momentumProblem, massProblem),
            std::make_tuple(momentumGridGeometry, massGridGeometry),
            std::make_tuple(momentumGridVariables, massGridVariables),
            couplingManager, timeLoop, xOld
        );

    // the linear solver
    using LinearSolver = LINEARSOLVER;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // the discrete L2 and Linfity errors
    const bool printErrors = getParam<bool>("Problem.PrintErrors", false);
    const bool printConvergenceTestFile = getParam<bool>("Problem.PrintConvergenceTestFile", false);
    NavierStokesTest::Errors errors(momentumProblem, massProblem, x);
    NavierStokesTest::ErrorCSVWriter errorCSVWriter(momentumProblem, massProblem);

    if (isStationary)
    {
        // linearize & solve
        Dune::Timer timer;
        nonLinearSolver.solve(x);

        // print discrete L2 and Linfity errors
        if (printErrors || printConvergenceTestFile)
        {
            errors.update(x);
            errorCSVWriter.printErrors(errors);

            if (printConvergenceTestFile)
                convergenceTestAppendErrors(momentumProblem, massProblem, errors);
        }

        // write vtk output
        analyticalSolVectors.update();
        vtkWriter.write(1.0);

        timer.stop();

        const auto& comm = Dune::MPIHelper::getCommunication();
        std::cout << "Simulation took " << timer.elapsed() << " seconds on "
                << comm.size() << " processes.\n"
                << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";
    }
    else
    {
        // time loop
        timeLoop->start(); do
        {
            // set the correct time level for the problem's boundary conditions
            momentumProblem->updateTime(timeLoop->time() + timeLoop->timeStepSize());
            massProblem->updateTime(timeLoop->time() + timeLoop->timeStepSize());

            // solve the non-linear system with time step control
            nonLinearSolver.solve(x, *timeLoop);

            // make the new solution the old solution
            xOld = x;
            momentumGridVariables->advanceTimeStep();
            massGridVariables->advanceTimeStep();

            // print discrete L2 and Linfity errors
            if (printErrors)
            {
                errors.update(x, timeLoop->time() + timeLoop->timeStepSize());
                errorCSVWriter.printErrors(errors);
            }

            // advance to the time loop to the next step
            timeLoop->advanceTimeStep();
            analyticalSolVectors.update(timeLoop->time());

            // write vtk output
            vtkWriter.write(timeLoop->time());

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
