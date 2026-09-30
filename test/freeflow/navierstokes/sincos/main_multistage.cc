// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the multi-stage multi-domain assembler with the coupled instationary CVFE
 *        Navier-Stokes mass and momentum models.
 *
 * The test checks that implicit Euler with the multi-stage multi-domain assembler reproduces
 * implicit Euler with the standard multi-domain assembler, that a single linear solve per stage
 * reproduces Newton's method for this linear problem, and that velocity and pressure
 * converge in time with the order of the scheme (self-convergence on a fixed grid).
 */

#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/pdesolver.hh>

#include <dumux/multidomain/assembler.hh>
#include <dumux/multidomain/multistagemultidomainassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

#include "properties.hh"

namespace Dumux::MultiStageTest {

using MomentumTypeTag = Properties::TTag::TYPETAG_MOMENTUM;
using MassTypeTag = Properties::TTag::TYPETAG_MASS;
using Scalar = GetPropType<MomentumTypeTag, Properties::Scalar>;
using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
using CouplingManager = GetPropType<MomentumTypeTag, Properties::CouplingManager>;
using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
using SolutionVector = typename Traits::SolutionVector;

constexpr auto momentumIdx = CouplingManager::freeFlowMomentumIndex;
constexpr auto massIdx = CouplingManager::freeFlowMassIndex;

struct Setup
{
    std::shared_ptr<CouplingManager> couplingManager;
    std::shared_ptr<MomentumProblem> momentumProblem;
    std::shared_ptr<MassProblem> massProblem;
    std::shared_ptr<MomentumGridVariables> momentumGridVariables;
    std::shared_ptr<MassGridVariables> massGridVariables;
    SolutionVector x;
    SolutionVector xOld;
};

std::unique_ptr<Setup> makeSetup(std::shared_ptr<const MomentumGridGeometry> momentumGridGeometry,
                                 std::shared_ptr<const MassGridGeometry> massGridGeometry)
{
    auto s = std::make_unique<Setup>();
    s->couplingManager = std::make_shared<CouplingManager>();
    s->momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, s->couplingManager);
    s->massProblem = std::make_shared<MassProblem>(massGridGeometry, s->couplingManager);
    s->momentumProblem->applyInitialSolution(s->x[momentumIdx]);
    s->massProblem->applyInitialSolution(s->x[massIdx]);
    s->xOld = s->x;

    s->momentumGridVariables = std::make_shared<MomentumGridVariables>(s->momentumProblem, momentumGridGeometry);
    s->massGridVariables = std::make_shared<MassGridVariables>(s->massProblem, massGridGeometry);
    s->couplingManager->init(s->momentumProblem, s->massProblem,
                             std::make_tuple(s->momentumGridVariables, s->massGridVariables), s->x, s->xOld);
    s->massGridVariables->init(s->x[massIdx]);
    s->momentumGridVariables->init(s->x[momentumIdx]);
    return s;
}

//! Implicit Euler with the standard multi-domain assembler and a fixed time step size
SolutionVector solveWithAssembler(std::shared_ptr<const MomentumGridGeometry> momentumGridGeometry,
                                  std::shared_ptr<const MassGridGeometry> massGridGeometry,
                                  Scalar dt, Scalar tEnd)
{
    auto s = makeSetup(momentumGridGeometry, massGridGeometry);
    auto timeLoop = std::make_shared<TimeLoop<Scalar>>(0.0, dt, tEnd, false);

    using Assembler = Experimental::MultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(s->momentumProblem, s->massProblem),
        std::make_tuple(momentumGridGeometry, massGridGeometry),
        std::make_tuple(s->momentumGridVariables, s->massGridVariables),
        s->couplingManager, timeLoop, s->xOld
    );

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager> nonLinearSolver(assembler, linearSolver, s->couplingManager);

    timeLoop->start(); do
    {
        s->momentumProblem->updateTime(timeLoop->time() + timeLoop->timeStepSize());
        s->massProblem->updateTime(timeLoop->time() + timeLoop->timeStepSize());
        nonLinearSolver.solve(s->x);
        s->xOld = s->x;
        s->momentumGridVariables->advanceTimeStep();
        s->massGridVariables->advanceTimeStep();
        timeLoop->advanceTimeStep();
    } while (!timeLoop->finished());

    return s->x;
}

//! Time integration with the multi-stage multi-domain assembler and a fixed time step size
template<bool singleLinearSolve = false>
SolutionVector solveWithMultiStageAssembler(std::shared_ptr<const MomentumGridGeometry> momentumGridGeometry,
                                            std::shared_ptr<const MassGridGeometry> massGridGeometry,
                                            std::shared_ptr<const Experimental::MultiStageMethod<Scalar>> method,
                                            Scalar dt, Scalar tEnd)
{
    auto s = makeSetup(momentumGridGeometry, massGridGeometry);

    using Assembler = Experimental::MultiStageMultiDomainAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(s->momentumProblem, s->massProblem),
        std::make_tuple(momentumGridGeometry, massGridGeometry),
        std::make_tuple(s->momentumGridVariables, s->massGridVariables),
        s->couplingManager, method, s->xOld
    );

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<Assembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    auto pdeSolver = [&]
    {
        if constexpr (singleLinearSolve)
            return std::make_shared<LinearPDESolver<Assembler, LinearSolver>>(assembler, linearSolver);
        else
            return std::make_shared<MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>>(
                assembler, linearSolver, s->couplingManager
            );
    }();

    using PDESolver = typename decltype(pdeSolver)::element_type;
    Experimental::MultiStageTimeStepper<PDESolver> timeStepper(pdeSolver, method);

    const auto numSteps = static_cast<int>(std::round(tEnd/dt));
    for (int stepIdx = 0; stepIdx < numSteps; ++stepIdx)
    {
        timeStepper.step(s->x, stepIdx*dt, dt);
        s->xOld = s->x;
        s->momentumGridVariables->advanceTimeStep();
        s->massGridVariables->advanceTimeStep();
    }

    return s->x;
}

template<class Vector>
Scalar discreteL2Norm(const Vector& x)
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

template<class Vector>
Scalar discreteL2Difference(const Vector& a, const Vector& b)
{
    auto diff = a;
    diff -= b;
    return discreteL2Norm(diff);
}

} // end namespace Dumux::MultiStageTest

int main(int argc, char** argv)
{
    using namespace Dumux;
    using namespace Dumux::MultiStageTest;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    GridManager<GetPropType<MomentumTypeTag, Properties::Grid>> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);

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
        const auto reference = solveWithAssembler(momentumGridGeometry, massGridGeometry, dt, tEnd);
        const auto multiStage = solveWithMultiStageAssembler(
            momentumGridGeometry, massGridGeometry,
            std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>(), dt, tEnd
        );
        const auto relDiffVelocity = discreteL2Difference(reference[momentumIdx], multiStage[momentumIdx])/discreteL2Norm(reference[momentumIdx]);
        const auto relDiffPressure = discreteL2Difference(reference[massIdx], multiStage[massIdx])/discreteL2Norm(reference[massIdx]);
        std::cout << "[Equivalence] implicit Euler: relative difference to the standard assembler: velocity "
                  << relDiffVelocity << ", pressure " << relDiffPressure
                  << " (tolerance " << equivalenceTolerance << ")" << std::endl;
        if (!(relDiffVelocity < equivalenceTolerance && relDiffPressure < equivalenceTolerance))
            passed = false;
    }

    using Method = Experimental::MultiStageMethod<Scalar>;

    // The problem is linear, so a single linear solve per stage has to reproduce Newton's method,
    // which requires the previous stages to enter the stage residual at their solutions. With a
    // numerically differentiated Jacobian, a single solve is exact only up to the finite-difference
    // error of the Jacobian, hence the separate tolerance.
    {
        const std::vector<std::shared_ptr<const Method>> multiStageMethods = {
            std::make_shared<Experimental::MultiStage::Theta<Scalar>>(0.5),
            std::make_shared<Experimental::MultiStage::DIRKSecondOrderAlexander<Scalar>>(),
            std::make_shared<Experimental::MultiStage::DIRKThirdOrderAlexander<Scalar>>()
        };

        const Scalar dt = tEnd/numStepsCoarse;
        for (const auto& method : multiStageMethods)
        {
            const auto newton = solveWithMultiStageAssembler(momentumGridGeometry, massGridGeometry, method, dt, tEnd);
            const auto linear = solveWithMultiStageAssembler<true>(momentumGridGeometry, massGridGeometry, method, dt, tEnd);
            const auto relDiffVelocity = discreteL2Difference(newton[momentumIdx], linear[momentumIdx])/discreteL2Norm(newton[momentumIdx]);
            const auto relDiffPressure = discreteL2Difference(newton[massIdx], linear[massIdx])/discreteL2Norm(newton[massIdx]);
            std::cout << "[Equivalence] " << method->name() << ": relative difference between a single linear solve"
                      << " per stage and Newton's method: velocity " << relDiffVelocity << ", pressure " << relDiffPressure
                      << " (tolerance " << singleSolveTolerance << ")" << std::endl;
            if (!(relDiffVelocity < singleSolveTolerance && relDiffPressure < singleSolveTolerance))
                passed = false;
        }
    }

    // temporal self-convergence of velocity and pressure
    const std::vector<std::tuple<std::shared_ptr<const Method>, int>> methods = {
        {std::make_shared<Experimental::MultiStage::ImplicitEuler<Scalar>>(), 1},
        {std::make_shared<Experimental::MultiStage::Theta<Scalar>>(0.5), 2}
    };

    for (const auto& [method, order] : methods)
    {
        std::vector<SolutionVector> results;
        for (int refIdx = 0; refIdx <= numRefinements; ++refIdx)
        {
            const Scalar dt = tEnd/(numStepsCoarse*(1 << refIdx));
            results.push_back(solveWithMultiStageAssembler(momentumGridGeometry, massGridGeometry, method, dt, tEnd));
        }

        std::vector<Scalar> velocityDifferences, pressureDifferences;
        for (std::size_t i = 0; i + 1 < results.size(); ++i)
        {
            velocityDifferences.push_back(discreteL2Difference(results[i][momentumIdx], results[i+1][momentumIdx]));
            pressureDifferences.push_back(discreteL2Difference(results[i][massIdx], results[i+1][massIdx]));
        }

        Scalar lastVelocityRate = 0.0, lastPressureRate = 0.0;
        for (std::size_t i = 0; i + 1 < velocityDifferences.size(); ++i)
        {
            lastVelocityRate = std::log2(velocityDifferences[i]/velocityDifferences[i+1]);
            lastPressureRate = std::log2(pressureDifferences[i]/pressureDifferences[i+1]);
            std::cout << "[Convergence] " << method->name() << ": velocity rate = " << lastVelocityRate
                      << ", pressure rate = " << lastPressureRate << std::endl;
        }

        const bool orderReached = lastVelocityRate > order - rateTolerance && lastPressureRate > order - rateTolerance;
        std::cout << "[Convergence] " << method->name() << ": expected order " << order
                  << (orderReached ? " reached" : " NOT reached") << std::endl;
        passed = passed && orderReached;
    }

    std::cout << (passed ? "Multi-stage multi-domain test: PASSED" : "Multi-stage multi-domain test: FAILED") << std::endl;
    return passed ? 0 : 1;
}
