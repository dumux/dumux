// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief Test of the backward differentiation formulas (BDF) and the multi-step time stepper
 *
 * - the variable-step weights reproduce the constant-step BDF tables, sum to zero and
 *   differentiate polynomials up to the order exactly,
 * - the number of time levels used for given step sizes follows the documented rules,
 * - the observed order of convergence for du/dt = -u + exp(t), u(0) = 1, with exact solution
 *   u(t) = cosh(t), for constant and variable step sizes (u''(0) is not zero, so that the order
 *   of the first steps shows in the global error),
 * - a step retried after a solver failure gives the same result as taking the reduced step directly,
 * - the constraint of the system u0' = -u0 + u1, 0 = sin(t) - u1 holds after every step.
 */
#include <config.h>

#include <cmath>
#include <array>
#include <vector>
#include <string>
#include <memory>
#include <random>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <functional>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dumux/io/format.hh>
#include <dumux/io/json.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/experimental/ode/odesolver.hh>
#include <dumux/experimental/timestepping/multistepmethods.hh>
#include <dumux/experimental/timestepping/multisteptimestepper.hh>

#include "scalarassembler.hh"
#include "constrainedode.hh"

namespace Dumux {

using Scalar = double;

void expect(const bool condition, const std::string& message)
{
    if (!condition)
        DUNE_THROW(Dune::InvalidStateException, message);
}

void expectNear(const Scalar value, const Scalar reference, const Scalar tolerance, const std::string& message)
{
    using std::abs;
    if (!(abs(value - reference) <= tolerance))
        DUNE_THROW(Dune::InvalidStateException, message << " (" << value << " vs. " << reference << ")");
}

//! Newton's method that reports a failure on the given call of apply, after modifying the variables
template<class Assembler, class LinearSolver>
class FailingNewton : public NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>
{
    using ParentType = NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>;
public:
    using ParentType::ParentType;
    using typename ParentType::Variables;

    void failOnCall(const std::size_t call)
    { failOnCall_ = call; }

    bool apply(Variables& vars) override
    {
        const bool converged = ParentType::apply(vars);
        return (++calls_ == failOnCall_) ? false : converged;
    }

private:
    std::size_t calls_ = 0;
    std::size_t failOnCall_ = 0;
};

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;
    Dumux::initialize(argc, argv);

    using BDF = Experimental::MultiStep::BDF<Scalar>;
    using MultiStepMethod = Experimental::MultiStepMethod<Scalar>;
    using MultiStageMethod = Experimental::MultiStageMethod<Scalar>;

    // constant-step BDF weights alpha_0, ..., alpha_k
    const std::vector<std::vector<Scalar>> tables = {
        {1.0, -1.0},
        {3.0/2.0, -2.0, 1.0/2.0},
        {11.0/6.0, -3.0, 3.0/2.0, -1.0/3.0},
        {25.0/12.0, -4.0, 3.0, -4.0/3.0, 1.0/4.0},
        {137.0/60.0, -5.0, 5.0, -10.0/3.0, 5.0/4.0, -1.0/5.0}
    };

    ////////////////////////////////////////////////////////////
    // weights
    ////////////////////////////////////////////////////////////
    for (std::size_t k = 1; k <= 5; ++k)
    {
        const BDF bdf(k);
        std::vector<Scalar> alpha(k+1), beta(k+1), stepSizes(k, 0.1);

        bdf.weights(stepSizes, alpha, beta);
        for (std::size_t j = 0; j <= k; ++j)
        {
            expectNear(alpha[j], tables[k-1][j], 1e-12, Fmt::format("Constant-step weight {} of BDF{}", j, k));
            expectNear(beta[j], j == 0 ? 1.0 : 0.0, 0.0, Fmt::format("Spatial weight {} of BDF{}", j, k));
        }

        std::mt19937 generator(k);
        std::uniform_real_distribution<Scalar> distribution(0.5, 1.5);
        for (int trial = 0; trial < 10; ++trial)
        {
            for (auto& dt : stepSizes)
                dt = 0.1*distribution(generator);
            bdf.weights(stepSizes, alpha, beta);

            // sum_j alpha_j p(tau_j) = dt p'(0) for p(tau) = tau^d, d <= k, with tau_0 = 0 the new time level
            std::vector<Scalar> tau(k+1, 0.0);
            for (std::size_t j = 1; j <= k; ++j)
                tau[j] = tau[j-1] - stepSizes[j-1];
            for (std::size_t d = 0; d <= k; ++d)
            {
                Scalar sum = 0.0, scale = stepSizes[0];
                for (std::size_t j = 0; j <= k; ++j)
                {
                    sum += alpha[j]*std::pow(tau[j], d);
                    scale += std::abs(alpha[j]*std::pow(tau[j], d));
                }
                expectNear(sum, d == 1 ? stepSizes[0] : 0.0, 1e-12*scale,
                           Fmt::format("BDF{} does not differentiate tau^{} exactly", k, d));
            }
        }
    }

    // closed form of the variable-step BDF2 weights with the step size ratio omega
    for (const Scalar omega : {0.25, 0.5, 1.0, 2.0})
    {
        const BDF bdf(2);
        std::vector<Scalar> alpha(3), beta(3);
        bdf.weights(std::vector<Scalar>{0.1*omega, 0.1}, alpha, beta);
        expectNear(alpha[0], (1.0 + 2.0*omega)/(1.0 + omega), 1e-13, "Variable-step BDF2 weight 0");
        expectNear(alpha[1], -(1.0 + omega), 1e-13, "Variable-step BDF2 weight 1");
        expectNear(alpha[2], omega*omega/(1.0 + omega), 1e-13, "Variable-step BDF2 weight 2");
    }

    ////////////////////////////////////////////////////////////
    // number of time levels used for given step sizes
    ////////////////////////////////////////////////////////////
    const auto numSteps = [](const MultiStepMethod& method, const std::vector<Scalar>& stepSizes)
    { return method.numStepsToUse(stepSizes); };

    expect(numSteps(BDF(5), {1.0, 1.0, 1.0, 1.0, 1.0}) == 5, "BDF5 with equal step sizes");
    expect(numSteps(BDF(5), {1.0, 1.0, 1.0}) == 3, "BDF5 with three available time levels");
    expect(numSteps(BDF(5), {1.0, 1.0, 1.5, 1.0, 1.0}) == 2, "BDF5 with a different third step size");
    expect(numSteps(BDF(5), {1.2, 1.0, 1.0, 1.0, 1.0}) == 2, "BDF5 with a larger next step");
    expect(numSteps(BDF(5), {3.0, 1.0, 1.0, 1.0, 1.0}) == 1, "BDF5 with a step size ratio above the BDF2 limit");
    expect(numSteps(BDF(2), {2.0, 1.0}) == 2, "BDF2 at the step size ratio limit");
    expect(numSteps(BDF(2), {2.0*(1.0 + 1e-12), 1.0}) == 1, "BDF2 above the step size ratio limit");
    expect(numSteps(BDF(2, 1.5), {1.6, 1.0}) == 1, "BDF2 with a smaller step size ratio limit");
    expect(numSteps(BDF(3), {1.0, 1.0, 1.0 + 1e-10}) == 3, "BDF3 with step sizes equal within the tolerance");

    bool thrown = false;
    try { BDF(2, 2.5); } catch (const Dune::InvalidStateException&) { thrown = true; }
    expect(thrown, "A step size ratio limit above 1 + sqrt(2) has to be rejected");

    ////////////////////////////////////////////////////////////
    // convergence for du/dt = -u + exp(t), u(0) = 1
    ////////////////////////////////////////////////////////////
    using Assembler = ScalarAssembler;
    using LinearSolver = ScalarLinearSolver;
    using Newton = NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>;
    using Variables = typename Assembler::Variables;

    const auto makeNewton = [](auto assembler)
    {
        return std::make_shared<Newton>(assembler, std::make_shared<LinearSolver>(),
                                        Dune::MPIHelper::getCommunication(), "", "Newton", 0);
    };

    const auto exact = [](const Scalar t) { using std::cosh; return cosh(t); };

    //! integrate with the given step sizes, returns the error at the end and the time levels used per step
    const auto integrate = [&](std::shared_ptr<const MultiStepMethod> method,
                               std::shared_ptr<const MultiStageMethod> startup,
                               const bool defaultStartup,
                               const std::vector<Scalar>& stepSizes)
    {
        auto newton = makeNewton(std::make_shared<Assembler>(-1.0));
        auto stepper = defaultStartup
            ? Experimental::MultiStepTimeStepper<Newton>(newton, method)
            : Experimental::MultiStepTimeStepper<Newton>(newton, method, startup);

        Variables vars(1.0);
        Scalar t = 0.0;
        std::vector<std::size_t> levels;
        for (const auto dt : stepSizes)
        {
            stepper.step(vars, t, dt);
            levels.push_back(stepper.lastNumSteps());
            t += dt;
        }

        using std::abs;
        return std::make_pair(abs(vars.dofs() - exact(t)), levels);
    };

    const auto constantSteps = [](const int n) { return std::vector<Scalar>(n, 1.0/n); };

    // smooth variation of the step size: t_i = g(i/n) with g(s) = s + a sin(2 pi s)/(2 pi)
    const auto variableSteps = [](const int n)
    {
        const Scalar a = 0.5;
        const auto g = [&](const Scalar s) { return s + a*std::sin(2.0*M_PI*s)/(2.0*M_PI); };
        std::vector<Scalar> steps(n);
        for (int i = 0; i < n; ++i)
            steps[i] = g(Scalar(i+1)/n) - g(Scalar(i)/n);
        return steps;
    };

    // constant steps with one step shortened to a tenth in the middle
    const auto shortenedStep = [](const int n)
    {
        std::vector<Scalar> steps(n, 1.0/n);
        steps.insert(steps.begin() + n/2, 0.1/n);
        return steps;
    };

    const std::vector<int> numStepsList = {10, 20, 40, 80, 160};

    const auto observedOrders = [&](auto&& errorForN)
    {
        std::vector<std::array<Scalar, 2>> convergence;
        for (const int n : numStepsList)
            convergence.push_back({1.0/n, errorForN(n)});
        std::vector<Scalar> orders;
        for (std::size_t i = 0; i + 1 < convergence.size(); ++i)
            orders.push_back(std::log2(convergence[i][1]/convergence[i+1][1]));
        return std::make_pair(convergence, orders);
    };

    const auto printOrders = [](const std::string& label, const std::vector<Scalar>& orders)
    {
        std::cout << "-- Observed orders " << label << ":";
        for (const auto o : orders)
            std::cout << Fmt::format(" {:.3f}", o);
        std::cout << std::endl;
    };

    Dumux::Json::JsonTree out;
    out["constant"] = Dumux::Json::JsonTree::object();
    out["variable"] = Dumux::Json::JsonTree::object();

    for (std::size_t k = 1; k <= 5; ++k)
    {
        const auto method = std::make_shared<BDF>(k);
        const auto [convergence, orders] = observedOrders([&](int n){ return integrate(method, nullptr, true, constantSteps(n)).first; });
        printOrders(Fmt::format("BDF{} with constant steps", k), orders);

        // the error of BDF5 approaches round-off for the smallest step sizes
        const auto checked = (k == 5) ? orders[orders.size()-2] : orders.back();
        expect(checked > k - 0.2, Fmt::format("Observed order of BDF{} with constant steps", k));

        auto& node = out["constant"][Fmt::format("bdf{}", k)];
        node["name"] = method->name();
        node["order"] = k;
        node["convergence"] = convergence;
    }

    // without a startup method of sufficient order, BDF3 reaches second order only
    {
        const auto method = std::make_shared<BDF>(3);
        const auto [convergence, orders] = observedOrders([&](int n){ return integrate(method, nullptr, false, constantSteps(n)).first; });
        printOrders("BDF3 started with BDF1 and BDF2", orders);
        expect(orders.back() > 1.8 && orders.back() < 2.3, "BDF3 started with lower BDF orders should be second order");

        auto& node = out["constant"]["bdf3_lowerorderstartup"];
        node["name"] = "BDF3 (started with BDF1, BDF2)";
        node["order"] = 2;
        node["convergence"] = convergence;
    }

    // variable step sizes: BDF1 and BDF2 keep their order, BDF3 falls back to BDF2
    for (std::size_t k = 1; k <= 3; ++k)
    {
        const auto method = std::make_shared<BDF>(k);
        const auto [convergence, orders] = observedOrders([&](int n){ return integrate(method, nullptr, true, variableSteps(n)).first; });
        printOrders(Fmt::format("BDF{} with variable steps", k), orders);
        const auto expectedOrder = std::min<std::size_t>(k, 2);
        expect(orders.back() > expectedOrder - 0.2 && orders.back() < expectedOrder + 0.3,
               Fmt::format("Observed order of BDF{} with variable steps", k));

        auto& node = out["variable"][Fmt::format("bdf{}", k)];
        node["name"] = method->name();
        node["order"] = expectedOrder;
        node["convergence"] = convergence;
    }

    // a step size ratio of ten after a shortened step: BDF2 takes that step with BDF1 and stays second order
    {
        const auto method = std::make_shared<BDF>(2);
        const auto [convergence, orders] = observedOrders([&](int n){ return integrate(method, nullptr, true, shortenedStep(n)).first; });
        printOrders("BDF2 with a shortened step", orders);
        expect(orders.back() > 1.8, "BDF2 with a shortened step should stay second order");

        const auto levels = integrate(method, nullptr, true, shortenedStep(20)).second;
        expect(levels[10] == 2 && levels[11] == 1 && levels[12] == 2,
               "BDF2 should take the step after the shortened step with one time level only");

        auto& node = out["variable"]["bdf2_shortenedstep"];
        node["name"] = "BDF2 (one step shortened to 1/10)";
        node["order"] = 2;
        node["convergence"] = convergence;
    }

    // the startup method is used while fewer time levels are available than the method uses
    {
        const auto levels = integrate(std::make_shared<BDF>(4), nullptr, true, constantSteps(10)).second;
        expect(levels[0] == 0 && levels[1] == 0 && levels[2] == 0 && levels[3] == 4,
               "BDF4 should take the first three steps with the startup method");
    }

    ////////////////////////////////////////////////////////////
    // a retried step gives the same result as the reduced step
    ////////////////////////////////////////////////////////////
    for (const std::size_t k : {2, 3})
    {
        const auto method = std::make_shared<BDF>(k);

        using Failing = FailingNewton<Assembler, LinearSolver>;
        auto failing = std::make_shared<Failing>(std::make_shared<Assembler>(-1.0), std::make_shared<LinearSolver>(),
                                                 Dune::MPIHelper::getCommunication(), "", "Newton", 0);
        // fail after the startup steps (two steps with five stages each for BDF3)
        failing->failOnCall(k == 2 ? 4 : 13);

        Experimental::MultiStepTimeStepper<Failing> stepper(failing, method);
        TimeLoop<Scalar> timeLoop(0.0, 0.1, 1.0, false);
        Variables vars(1.0);
        std::vector<Scalar> stepSizes;
        while (!timeLoop.finished())
        {
            stepper.step(vars, timeLoop);
            stepSizes.push_back(timeLoop.timeStepSize());
            timeLoop.advanceTimeStep();
        }

        expect(std::count_if(stepSizes.begin(), stepSizes.end(), [](Scalar dt){ return dt < 0.075; }) > 0,
               "The failure should have reduced the step size");

        auto newton = makeNewton(std::make_shared<Assembler>(-1.0));
        Experimental::MultiStepTimeStepper<Newton> reference(newton, method);
        Variables refVars(1.0);
        Scalar t = 0.0;
        for (const auto dt : stepSizes)
        {
            reference.step(refVars, t, dt);
            t += dt;
        }

        expectNear(vars.dofs(), refVars.dofs(), 1e-14,
                   Fmt::format("BDF{} with a retried step differs from taking the reduced step directly", k));
        std::cout << Fmt::format("-- BDF{} with a retried step: difference to the reduced step {:.2e}\n",
                                 k, std::abs(vars.dofs() - refVars.dofs()));
    }

    ////////////////////////////////////////////////////////////
    // the constraint u1 = sin(t) holds after every step
    ////////////////////////////////////////////////////////////
    {
        using ODEAssembler = Experimental::MultiStageODEAssembler<Test::ConstrainedODE>;
        using ODELinearSolver = Experimental::ODELocalLinearSolver<Test::ConstrainedODE::JacobianMatrix, Test::ConstrainedODE::ResidualType>;
        using ODENewton = NewtonSolver<ODEAssembler, ODELinearSolver, DefaultPartialReassembler>;
        using ODEVars = Test::ConstrainedODE::Variables;

        for (std::size_t k = 1; k <= 5; ++k)
        {
            const auto method = std::make_shared<BDF>(k);
            Scalar maxConstraintError = 0.0;

            const auto [convergence, orders] = observedOrders([&](int n)
            {
                auto newton = std::make_shared<ODENewton>(std::make_shared<ODEAssembler>(std::make_shared<const Test::ConstrainedODE>()),
                                                          std::make_shared<ODELinearSolver>(),
                                                          Dune::MPIHelper::getCommunication(), "", "Newton", 0);
                Experimental::MultiStepTimeStepper<ODENewton> stepper(newton, method);
                ODEVars vars(Test::ConstrainedODE::SolutionVector{0.0, 0.0});
                Scalar t = 0.0;
                for (const auto dt : constantSteps(n))
                {
                    stepper.step(vars, t, dt);
                    t += dt;
                    using std::abs; using std::max; using std::sin;
                    maxConstraintError = max(maxConstraintError, abs(vars.dofs()[1] - sin(t)));
                }
                using std::abs;
                return abs(vars.dofs()[0] - Test::constrainedExact(t));
            });

            printOrders(Fmt::format("BDF{} for the constrained system", k), orders);
            std::cout << Fmt::format("-- Largest constraint violation with BDF{}: {:.2e}\n", k, maxConstraintError);
            expect(maxConstraintError < 1e-13, Fmt::format("BDF{} violates the constraint", k));
            const auto checked = (k == 5) ? orders[orders.size()-2] : orders.back();
            expect(checked > k - 0.2, Fmt::format("Observed order of BDF{} for the constrained system", k));
        }
    }

    const std::string outputFile = "test_multistepmethods_data.json";
    std::ofstream output(outputFile);
    output << std::setw(2) << out << std::endl;
    std::cout << "Wrote convergence data to " << outputFile << std::endl;

    return 0;
}
