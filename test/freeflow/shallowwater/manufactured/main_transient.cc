// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Time-dependent manufactured solution of the shallow water equations integrated with a BDF method:
 *        errors with respect to the exact solution at the final time
 *
 * The BDF order is `TimeStepper.BDFOrder` (default 2). The time step size is constant, `TimeLoop.Dt` rounded
 * such that an integer number of steps reaches `TimeLoop.TEnd`.
 */
#include <config.h>

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <type_traits>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#if HAVE_DUNE_UGGRID
#include <dumux/io/grid/gridmanager_ug.hh>
#endif
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/experimental/assembly/multistagefvassembler.hh>
#include <dumux/experimental/timestepping/multistepmethods.hh>
#include <dumux/experimental/timestepping/multisteptimestepper.hh>

#include <test/experimental/timestepping/implicitsolvecheck.hh>
#include <test/freeflow/shallowwater/errors.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();
    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::TYPETAG;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto numSteps = static_cast<std::size_t>(std::round(tEnd/getParam<Scalar>("TimeLoop.Dt")));
    const Scalar dt = tEnd/numSteps;

    const auto order = getParam<std::size_t>("TimeStepper.BDFOrder", 2);
    auto method = std::make_shared<Experimental::MultiStep::BDF<Scalar>>(order);

    using Assembler = Experimental::MultiStageFVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, method, xOld);

    using LinearSolver = LINEARSOLVER<LinearSolverTraits<GridGeometry>,
                                      LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = [&]{
        if constexpr (std::is_constructible_v<LinearSolver, decltype(leafGridView), decltype(gridGeometry->dofMapper())>)
            return std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());
        else
            return std::make_shared<LinearSolver>();
    }();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver);

    Experimental::MultiStepTimeStepper<NewtonSolver> timeStepper(nonLinearSolver, method);

    Test::checkImplicitSolveContract(*assembler, x, 0.0, 1.0/dt, 1.0);

    // number of steps by the number of previous time levels used (zero for a startup step)
    std::array<std::size_t, 6> numStepsByLevels{};
    for (std::size_t n = 0; n < numSteps; ++n)
    {
        timeStepper.step(x, n*dt, dt);
        ++numStepsByLevels[timeStepper.lastNumSteps()];
        xOld = x;
        gridVariables->advanceTimeStep();
    }

    if (mpiHelper.rank() == 0)
    {
        std::cout << numSteps << " steps of size " << dt << " with " << method->name() << ", by number of previous time levels:";
        for (std::size_t k = 0; k <= order; ++k)
            std::cout << " " << k << ": " << numStepsByLevels[k];
        std::cout << std::endl;
    }

    using Solution = ShallowWaterTest::ManufacturedSolution<Scalar>;
    const auto errors = ShallowWaterTest::computeErrors(*gridGeometry, x, [&](const auto& element){
        const auto [h, u, v, z] = Solution::state(element.geometry().center(), tEnd, problem->frequency());
        return std::array<Scalar, 3>{{ h, u, v }};
    });
    if (mpiHelper.rank() == 0)
        ShallowWaterTest::printErrors(errors, tEnd);

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
