// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the multi-stage assembler with the instationary CVFE momentum model.
 *
 * The momentum balance is solved with the analytical pressure of the instationary sincos
 * problem. The test checks that implicit Euler with the multi-stage assembler reproduces
 * implicit Euler with the standard assembler, that the residual assembled without the Jacobian
 * equals the one assembled with it, that a single linear solve per stage reproduces Newton's
 * method for this linear problem, and that the implicit multi-stage schemes converge in time
 * with their theoretical order (self-convergence on a fixed grid).
 */

#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/pdesolver.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/assembly/assembler.hh>
#include <dumux/assembly/multistageassembler.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

#include <test/freeflow/navierstokes/errors_cvfe.hh>

#include "properties.hh"

#ifndef TYPETAG_MOMENTUM_ONLY
#define TYPETAG_MOMENTUM_ONLY SincosTestMomentumOnlyPQ1BubbleHybrid
#endif

namespace Dumux {

using TypeTag = Properties::TTag::TYPETAG_MOMENTUM_ONLY;
using Scalar = GetPropType<TypeTag, Properties::Scalar>;
using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
using Problem = GetPropType<TypeTag, Properties::Problem>;
using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

struct Result
{
    SolutionVector x;
    std::shared_ptr<Problem> problem;
    std::shared_ptr<GridVariables> gridVariables;
};

//! Implicit Euler with the standard assembler and a fixed time step size
Result solveWithAssembler(std::shared_ptr<const GridGeometry> gridGeometry, Scalar dt, Scalar tEnd)
{
    auto problem = std::make_shared<Problem>(gridGeometry);
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd, false);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, timeLoop, xOld);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, linearSolver);

    timeLoop->start(); do
    {
        problem->updateTime(timeLoop->time() + timeLoop->timeStepSize());
        nonLinearSolver.solve(x);
        xOld = x;
        gridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();
    } while (!timeLoop->finished());

    return {x, problem, gridVariables};
}

//! Time integration with the multi-stage assembler and a fixed time step size
template<bool singleLinearSolve = false>
Result solveWithMultiStageAssembler(std::shared_ptr<const GridGeometry> gridGeometry,
                                    std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> method,
                                    Scalar dt, Scalar tEnd)
{
    auto problem = std::make_shared<Problem>(gridGeometry);
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = Experimental::MultiStageAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, method, xOld);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    using PDESolver = std::conditional_t<singleLinearSolve,
                                         LinearPDESolver<Assembler, LinearSolver>,
                                         NewtonSolver<Assembler, LinearSolver>>;
    auto pdeSolver = std::make_shared<PDESolver>(assembler, linearSolver);

    Experimental::MultiStageTimeStepper<PDESolver> timeStepper(pdeSolver, method);

    const auto numSteps = static_cast<int>(std::round(tEnd/dt));
    for (int stepIdx = 0; stepIdx < numSteps; ++stepIdx)
    {
        timeStepper.step(x, stepIdx*dt, dt);
        xOld = x;
        gridVariables->advanceTimeStep();
    }

    problem->setTime(tEnd);
    return {x, problem, gridVariables};
}

/*!
 * \brief Largest relative difference over the stages of one time step between the residual
 *        assembled without and with the Jacobian, both at the solution the stage starts from
 */
Scalar residualAssemblyDifference(std::shared_ptr<const GridGeometry> gridGeometry,
                                  std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> method,
                                  Scalar dt)
{
    auto problem = std::make_shared<Problem>(gridGeometry);
    SolutionVector x;
    problem->applyInitialSolution(x);
    auto xOld = x;

    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = Experimental::MultiStageAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables, method, xOld);

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    NewtonSolver<Assembler, LinearSolver> nonLinearSolver(assembler, std::make_shared<LinearSolver>());

    Scalar maxRelDiff = 0.0;
    for (std::size_t stageIdx = 1; stageIdx <= method->numStages(); ++stageIdx)
    {
        assembler->prepareStage(x, std::make_shared<Experimental::MultiStageParams<Scalar>>(*method, stageIdx, 0.0, dt));

        assembler->assembleJacobianAndResidual(x);
        auto difference = assembler->residual();
        assembler->assembleResidual(x);
        difference -= assembler->residual();

        using std::max;
        maxRelDiff = max(maxRelDiff, difference.two_norm()/assembler->residual().two_norm());

        nonLinearSolver.solve(x);
    }

    return maxRelDiff;
}

Scalar discreteL2Norm(const SolutionVector& x)
{
    Scalar sum = 0.0;
    std::size_t n = 0;
    for (const auto& block : x)
        for (const auto& value : block)
        {
            sum += value*value;
            ++n;
        }
    return std::sqrt(sum/n);
}

Scalar discreteL2Difference(const SolutionVector& a, const SolutionVector& b)
{
    auto diff = a;
    diff -= b;
    return discreteL2Norm(diff);
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    GridManager<GetPropType<TypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());

    const auto tEnd = getParam<Scalar>("TimeLoop.TEnd");
    const auto numStepsCoarse = getParam<int>("MultiStageTest.NumStepsCoarse");
    const auto numRefinements = getParam<int>("MultiStageTest.NumTimeStepRefinements");
    const auto equivalenceTolerance = getParam<Scalar>("MultiStageTest.EquivalenceTolerance");
    const auto singleSolveTolerance = getParam<Scalar>("MultiStageTest.SingleSolveTolerance");
    const auto rateTolerance = getParam<Scalar>("MultiStageTest.RateTolerance");

    bool passed = true;

    // implicit Euler with the multi-stage assembler has to reproduce the standard assembler
    {
        const Scalar dt = tEnd/numStepsCoarse;
        const auto reference = solveWithAssembler(gridGeometry, dt, tEnd);
        const auto multiStage = solveWithMultiStageAssembler(
            gridGeometry, std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>(), dt, tEnd
        );
        const auto relDiff = discreteL2Difference(reference.x, multiStage.x)/discreteL2Norm(reference.x);
        std::cout << "[Equivalence] implicit Euler: relative difference to the standard assembler = "
                  << relDiff << " (tolerance " << equivalenceTolerance << ")" << std::endl;
        if (!(relDiff < equivalenceTolerance))
            passed = false;
    }

    using Method = Experimental::MultiStageMethod<Scalar>;
    const std::vector<std::tuple<std::shared_ptr<const Method>, int>> methods = {
        {std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>(), 1},
        {std::make_shared<Experimental::MultiStage::Theta<Scalar>>(0.5), 2},
        {std::make_shared<Experimental::MultiStage::DIRKSecondOrderAlexander<Scalar>>(), 2},
        {std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<Scalar>>(), 3}
    };

    {
        const auto method = std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<Scalar>>();
        const auto relDiff = residualAssemblyDifference(gridGeometry, method, tEnd/numStepsCoarse);
        std::cout << "[Equivalence] " << method->name() << ": largest relative difference between the residual"
                  << " assembled without and with the Jacobian = " << relDiff << std::endl;
        if (!(relDiff < 1e-12))
            passed = false;
    }

    // The problem is linear, so a single linear solve per stage has to reproduce Newton's method,
    // which requires the previous stages to enter the stage residual at their solutions. With a
    // numerically differentiated Jacobian, a single solve is exact only up to the finite-difference
    // error of the Jacobian, hence the separate tolerance.
    for (const auto& [method, order] : methods)
    {
        const Scalar dt = tEnd/numStepsCoarse;
        const auto newton = solveWithMultiStageAssembler(gridGeometry, method, dt, tEnd);
        const auto linear = solveWithMultiStageAssembler<true>(gridGeometry, method, dt, tEnd);
        const auto relDiff = discreteL2Difference(newton.x, linear.x)/discreteL2Norm(newton.x);
        std::cout << "[Equivalence] " << method->name() << ": relative difference between a single linear solve"
                  << " per stage and Newton's method = " << relDiff << " (tolerance " << singleSolveTolerance << ")" << std::endl;
        if (!(relDiff < singleSolveTolerance))
            passed = false;
    }

    // temporal self-convergence of the implicit multi-stage schemes
    for (const auto& [method, order] : methods)
    {
        std::vector<Result> results;
        for (int refIdx = 0; refIdx <= numRefinements; ++refIdx)
        {
            const Scalar dt = tEnd/(numStepsCoarse*(1 << refIdx));
            results.push_back(solveWithMultiStageAssembler(gridGeometry, method, dt, tEnd));

            const auto [volume, errors] = calculateL2AndH1Errors(*results.back().problem, *results.back().gridVariables, results.back().x);
            std::cout << "[Convergence] " << method->name() << ": dt = " << dt
                      << " velocity L2 error to the analytical solution = " << errors[0] << std::endl;
        }

        std::vector<Scalar> differences;
        for (std::size_t i = 0; i + 1 < results.size(); ++i)
            differences.push_back(discreteL2Difference(results[i].x, results[i+1].x));

        Scalar lastRate = 0.0;
        for (std::size_t i = 0; i + 1 < differences.size(); ++i)
        {
            lastRate = std::log2(differences[i]/differences[i+1]);
            std::cout << "[Convergence] " << method->name() << ": rate = " << lastRate << std::endl;
        }

        const bool orderReached = lastRate > order - rateTolerance;
        std::cout << "[Convergence] " << method->name() << ": expected order " << order
                  << (orderReached ? " reached" : " NOT reached") << std::endl;
        passed = passed && orderReached;
    }

    std::cout << (passed ? "Multi-stage momentum test: PASSED" : "Multi-stage momentum test: FAILED") << std::endl;
    return passed ? 0 : 1;
}
