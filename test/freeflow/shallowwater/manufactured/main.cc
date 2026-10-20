// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterTests
 * \brief Steady manufactured solution of the shallow water equations: errors with respect to the exact solution
 *
 * The steady state is computed directly with the Newton method, or, if `TimeLoop.TEnd` is given, by implicit
 * Euler steps of growing size until that time.
 */
#include <config.h>

#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <type_traits>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#if HAVE_DUNE_UGGRID
#include <dumux/io/grid/gridmanager_ug.hh>
#endif
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/assembly/fvassembler.hh>

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

    // the sources agree with an independent symbolic derivation (sympy) at points (position, time, frequency)
    {
        using Solution = ShallowWaterTest::ManufacturedSolution<Scalar>;
        struct Reference { Dune::FieldVector<Scalar, 2> pos; Scalar t, omega; std::array<Scalar, 3> values; };
        const std::array<Reference, 4> reference{{
            {{0.3, 0.7}, 0.0, 0.0, {-1.294671721345872e-01, 2.375026512891299e+00, -5.207677065976531e+00}},
            {{0.81, 0.12}, 0.0, 0.0, {1.085241837853164e-01, -3.674026563380908e+00, 6.748216636764020e-01}},
            {{0.3, 0.7}, 0.4, 2.0, {2.700309032355022e-02, 3.787150641898493e+00, -5.740980113906270e+00}},
            {{0.81, 0.12}, 1.3, 2.0, {-8.603161464105410e-02, 1.985199101188140e-01, 2.226299324849193e+00}}
        }};
        for (const auto& [pos, t, omega, values] : reference)
        {
            const auto source = Solution::source(pos, 9.81, t, omega);
            for (int k = 0; k < 3; ++k)
                if (std::abs(source[k] - values[k]) > 1e-12*(1.0 + std::abs(values[k])))
                    DUNE_THROW(Dune::Exception, "Source " << k << " at " << pos << ", t = " << t << ", omega = " << omega
                                                << " is " << source[k] << ", expected " << values[k]);
        }
    }

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

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    const bool pseudoTransient = hasParam("TimeLoop.TEnd");
    auto xOld = x;
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, getParam<Scalar>("TimeLoop.DtInitial", 1.0),
                                                       getParam<Scalar>("TimeLoop.TEnd", 1.0));
    timeLoop->setMaxTimeStepSize(getParam<Scalar>("TimeLoop.MaxTimeStepSize", 1.0));
    auto assembler = pseudoTransient ? std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld)
                                     : std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using LinearSolver = LINEARSOLVER<LinearSolverTraits<GridGeometry>,
                                      LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = [&]{
        if constexpr (std::is_constructible_v<LinearSolver, decltype(leafGridView), decltype(gridGeometry->dofMapper())>)
            return std::make_shared<LinearSolver>(leafGridView, gridGeometry->dofMapper());
        else
            return std::make_shared<LinearSolver>();
    }();

    using NewtonSolver = Dumux::NewtonSolver<Assembler, LinearSolver>;
    NewtonSolver nonLinearSolver(assembler, linearSolver);
    if (pseudoTransient)
    {
        timeLoop->start(); do
        {
            nonLinearSolver.solve(x, *timeLoop);
            xOld = x;
            gridVariables->advanceTimeStep();
            timeLoop->advanceTimeStep();
            timeLoop->reportTimeStep();
            timeLoop->setTimeStepSize(nonLinearSolver.suggestTimeStepSize(timeLoop->timeStepSize()));
        } while (!timeLoop->finished());
    }
    else
        nonLinearSolver.solve(x);

    using Solution = ShallowWaterTest::ManufacturedSolution<Scalar>;
    const auto errors = ShallowWaterTest::computeErrors(*gridGeometry, x, [&](const auto& element){
        const auto [h, u, v, z] = Solution::state(element.geometry().center());
        return std::array<Scalar, 3>{{ h, u, v }};
    });
    if (mpiHelper.rank() == 0)
        ShallowWaterTest::printErrors(errors);

    if (getParam<bool>("Problem.EnableVtkOutput", false))
    {
        using IOFields = GetPropType<TypeTag, Properties::IOFields>;
        VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
        IOFields::initOutputModule(vtkWriter);
        vtkWriter.write(1.0);
    }

    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
