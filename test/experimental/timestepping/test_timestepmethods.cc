//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <cmath>
#include <cassert>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <memory>
#include <optional>
#include <array>
#include <string>
#include <vector>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>

#include <dumux/io/format.hh>
#include <dumux/io/json.hh>
#include <dumux/common/initialize.hh>
#include <dumux/experimental/common/variables.hh>
#include <dumux/nonlinear/newtonsolver.hh>

#include <dumux/experimental/timestepping/timelevel.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

/*
   This tests the time integration methods by solving the
   linear ODE: du/dt - exp(t) = 0, where u is the unknown
   and t is the time. We use the initial condition u_0 = 1,
   and thus, the exact solution is u_e = exp(t) - 1.
 */

namespace Dumux {

class ScalarAssembler
{
public:
    using Scalar = double;
    using SolutionVector = Scalar;
    using ResidualType = Scalar;
    using JacobianMatrix = Scalar;
    using Variables = Experimental::Variables<SolutionVector>;
    using SolveParams = Experimental::ImplicitSolveParams<Scalar, ResidualType>;

    explicit ScalarAssembler(Scalar lambda = 0.0) : lambda_(lambda) {}

    void setLinearSystem() {}
    JacobianMatrix& jacobian() { return jac_; }
    ResidualType& residual() { return res_; }

    void assembleResidual(const Variables& vars)
    {
        if (!params_)
            DUNE_THROW(Dune::InvalidStateException, "No implicit solve set");

        const auto& weights = params_->weights;
        res_ = weights.storage*storage_(vars) + weights.spatial*spatial_(vars);
        if (params_->constant)
            res_ += *params_->constant;
    }

    void assembleJacobianAndResidual(const Variables& vars)
    {
        assembleResidual(vars);
        jac_ = params_->weights.storage - lambda_*params_->weights.spatial;
    }

    void setImplicitSolve(Variables& vars, const SolveParams& params)
    {
        vars.updateTime(params.timeLevel);
        params_ = params;
    }

    void assembleOperators(Variables& vars,
                           const Experimental::TimeLevel<Scalar>& timeLevel,
                           Experimental::OperatorEvaluations<ResidualType> evaluations)
    {
        vars.updateTime(timeLevel);
        if (evaluations.storage)
            *evaluations.storage = storage_(vars);
        if (evaluations.spatial)
            *evaluations.spatial = spatial_(vars);
    }

private:
    Scalar storage_(const Variables& vars) const
    { return vars.dofs(); }

    Scalar spatial_(const Variables& vars) const
    {
        using std::exp;
        return -(lambda_*vars.dofs() + exp(vars.timeLevel().current()));
    }

    Scalar lambda_;
    ResidualType res_;
    JacobianMatrix jac_;
    std::optional<SolveParams> params_;
};

class ScalarLinearSolver
{
public:
    void setResidualReduction(double residualReduction) {}

    bool solve(const double& A, double& x, const double& b) const
    { x = b/A; return true; }

    double norm(const double residual) const
    {
        using std::abs;
        return abs(residual);
    }
};

} // end namespace Dumux

int main(int argc, char* argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    using Assembler = ScalarAssembler;
    using LinearSolver = ScalarLinearSolver;
    using NewtonSolver = NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler>;

    auto assembler = std::make_shared<Assembler>();
    auto linearSolver = std::make_shared<LinearSolver>();
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver);

    // initial solution and variables
    using Scalar = typename Assembler::Scalar;
    using Variables = typename Assembler::Variables;
    using SolutionVector = typename Variables::SolutionVector;

    const auto exact = [] (const Scalar t)
    {
        using std::exp;
        return exp(t) - 1.0;
    };

    const auto computeError = [&] (const Variables& vars)
    {
        using std::abs;
        const auto time = vars.timeLevel().current();
        const auto exactSol = exact(time);
        const auto absErr = abs(vars.dofs()-exactSol);
        const auto relErr = absErr/exactSol;
        return std::make_pair(absErr, relErr);
    };

    using namespace Experimental::MultiStage;

    const auto expectNear = [] (const auto value,
                                const auto reference,
                                const auto tolerance,
                                const std::string& message)
    {
        using std::abs;
        if (abs(value - reference) > tolerance)
            DUNE_THROW(Dune::InvalidStateException, message);
    };

    const auto integrateTo = [&] (const auto& method, const Scalar dt, const Scalar endTime = 1.0)
    {
        SolutionVector x = 0.0;
        Variables vars(x);

        using TimeStepper = Experimental::MultiStageTimeStepper<NewtonSolver>;
        TimeStepper timeStepper(newtonSolver, method);

        Scalar time = 0.0;
        while (Dune::FloatCmp::lt(time, endTime, 1e-12))
        {
            const auto stepDt = std::min(dt, endTime - time);
            timeStepper.step(vars, time, stepDt);
            time += stepDt;
        }

        return std::make_pair(vars, computeError(vars));
    };

    const auto testMethodProperties = [&] (const auto& method, const bool expectedImplicit)
    {
        using std::abs;
        constexpr Scalar eps = 1e-12;

        if (method->implicit() != expectedImplicit)
            DUNE_THROW(Dune::InvalidStateException,
                       "Unexpected implicit/explicit classification for " + method->name());

        expectNear(method->timeStepWeight(0), 0.0, eps,
                   "The first stage time weight has to be zero for " + method->name());
        expectNear(method->timeStepWeight(method->numStages()), 1.0, eps,
                   "The final stage has to advance to t^{n+1} for " + method->name());

        Scalar maxDiagonalWeight = 0.0;
        for (std::size_t i = 1; i <= method->numStages(); ++i)
        {
            Scalar alphaSum = 0.0;
            for (std::size_t k = 0; k <= i; ++k)
                alphaSum += method->temporalWeight(i, k);

            expectNear(alphaSum, 0.0, eps,
                       "Temporal weights do not preserve constant states for " + method->name());

            const auto diagonalWeight = method->spatialWeight(i, i);
            maxDiagonalWeight = std::max(maxDiagonalWeight, diagonalWeight);
            if (expectedImplicit)
            {
                // DIRK stages have a non-negative diagonal weight; a zero diagonal marks an
                // explicit stage, e.g. the appended (non stiffly-accurate) update stage of the
                // Qin-Zhang scheme
                if (diagonalWeight < 0.0)
                    DUNE_THROW(Dune::InvalidStateException,
                               "DIRK stage weights must be non-negative for " + method->name());
            }
            else
                expectNear(diagonalWeight, 0.0, eps,
                           "Explicit methods must have zero diagonal stage weights for " + method->name());
        }

        // an implicit method has to have at least one genuinely implicit (positive diagonal) stage
        if (expectedImplicit && !(maxDiagonalWeight > 0.0))
            DUNE_THROW(Dune::InvalidStateException,
                       "Implicit methods need at least one positive diagonal stage weight for " + method->name());
    };

    const auto testSingleStepAccuracy = [&] (const auto& method, const Scalar maxRelError)
    {
        std::cout << "\n-- Integration with " << method->name() << ":\n\n";
        const auto [vars, error] = integrateTo(method, 0.01, 0.01);
        const auto [abs, rel] = error;
        std::cout << "\n"
                  << "-- Summary\n"
                  << "-- ===========================\n"
                  << "-- Exact solution:    " << Fmt::format("{:.4e}\n", exact(0.01))
                  << "-- Discrete solution: " << Fmt::format("{:.4e}\n", vars.dofs())
                  << "-- Absolute error:    " << Fmt::format("{:.4e}\n", abs)
                  << "-- Relative error:    " << Fmt::format("{:.4e}", rel) << std::endl;

        if (!(rel < maxRelError))
            DUNE_THROW(Dune::InvalidStateException,
                       "Single-step error is larger than expected for " + method->name());
    };

    const auto testObservedOrder = [&] (const auto& method, const int expectedOrder)
    {
        using std::log2;

        const auto coarseError = integrateTo(method, 0.25).second.first;
        const auto fineError = integrateTo(method, 0.125).second.first;
        const auto observedOrder = log2(coarseError/fineError);

        std::cout << "-- Observed order for " << method->name() << ": "
                  << observedOrder << std::endl;

        if (observedOrder + 0.2 < expectedOrder)
            DUNE_THROW(Dune::InvalidStateException,
                       "Observed order is lower than expected for " + method->name());
    };

    using Method = std::shared_ptr<const Experimental::MultiStageMethod<Scalar>>;
    const Method explicitEuler = std::make_shared<ExplicitEuler<Scalar>>();
    const Method implicitEuler = std::make_shared<ImplicitEuler<Scalar>>();
    const Method crankNicolson = std::make_shared<Theta<Scalar>>(0.5);
    const Method heun = std::make_shared<RungeKuttaExplicitSecondOrderHeun<Scalar>>();
    const Method rk3 = std::make_shared<RungeKuttaExplicitThirdOrder<Scalar>>();
    const Method rk4 = std::make_shared<RungeKuttaExplicitFourthOrder<Scalar>>();
    const Method dirk2 = std::make_shared<DIRKSecondOrderAlexander<Scalar>>();
    const Method dirk3 = std::make_shared<DIRKThirdOrderAlexander<Scalar>>();
    const Method qinZhang = std::make_shared<QinZhangSymplecticDIRK<Scalar>>();

    // shorter display names for plot legends
    const auto displayName = [](const std::string& id, const Method& method)
    {
        if (id == "crank_nicolson") return std::string{"Crank-Nicolson"};
        if (id == "dirk2") return std::string{"DIRK2 (Alexander)"};
        if (id == "dirk3") return std::string{"DIRK3 (Alexander)"};
        if (id == "qin_zhang") return std::string{"Qin-Zhang (symplectic)"};
        return method->name();
    };

    const std::vector<std::pair<Method, bool>> methodsToCheck = {
        {explicitEuler, false},
        {implicitEuler, true},
        {crankNicolson, true},
        {heun, false},
        {rk3, false},
        {rk4, false},
        {dirk2, true},
        {dirk3, true},
        {qinZhang, true}
    };

    for (const auto& [method, expectedImplicit] : methodsToCheck)
        testMethodProperties(method, expectedImplicit);

    // when adding a new method here, run first to verify approximate error
    // and then add here as tolerance for the single-step accuracy regression test
    testSingleStepAccuracy(explicitEuler, 6e-3);
    testSingleStepAccuracy(implicitEuler, 6e-3);
    testSingleStepAccuracy(crankNicolson, 1e-5);
    testSingleStepAccuracy(heun, 1e-5);
    testSingleStepAccuracy(rk3, 3e-8);
    testSingleStepAccuracy(rk4, 1e-11);
    testSingleStepAccuracy(dirk2, 1e-5);
    testSingleStepAccuracy(dirk3, 1e-8);
    testSingleStepAccuracy(qinZhang, 2e-6);

    testObservedOrder(explicitEuler, 1);
    testObservedOrder(implicitEuler, 1);
    testObservedOrder(crankNicolson, 2);
    testObservedOrder(heun, 2);
    testObservedOrder(rk3, 3);
    testObservedOrder(rk4, 4);
    testObservedOrder(dirk2, 2);
    testObservedOrder(dirk3, 3);
    testObservedOrder(qinZhang, 2);

    // --- convergence sweep: error vs. step size for all methods, written to JSON for plotting ---
    // Use a second test equation, du/dt = -u + exp(t), u(0) = 0, with exact solution u(t) = sinh(t).
    // Unlike du/dt = exp(t) above, this right-hand side depends on u, so the methods no longer
    // collapse onto shared quadrature rules (cf. the Crank-Nicolson/Heun and RK3/RK4 degeneracy
    // for a right-hand side that depends on t only).
    const Scalar lambda = -1.0;
    auto assemblerLin = std::make_shared<Assembler>(lambda);
    auto newtonSolverLin = std::make_shared<NewtonSolver>(assemblerLin, linearSolver);

    const auto exactLin = [] (const Scalar t)
    {
        using std::sinh;
        return sinh(t);
    };

    const auto integrateToLin = [&] (const auto& method, const Scalar dt, const Scalar endTime = 1.0)
    {
        SolutionVector x = 0.0;
        Variables vars(x);

        Experimental::MultiStageTimeStepper<NewtonSolver> timeStepper(newtonSolverLin, method);

        Scalar time = 0.0;
        while (Dune::FloatCmp::lt(time, endTime, 1e-12))
        {
            const auto stepDt = std::min(dt, endTime - time);
            timeStepper.step(vars, time, stepDt);
            time += stepDt;
        }

        using std::abs;
        return abs(vars.dofs() - exactLin(vars.timeLevel().current()));
    };

    struct MethodInfo
    {
        std::string id;
        Method method;
        int order;
    };

    const std::vector<MethodInfo> allMethods = {
        {"explicit_euler", explicitEuler, 1},
        {"implicit_euler", implicitEuler, 1},
        {"crank_nicolson", crankNicolson, 2},
        {"heun",           heun,          2},
        {"dirk2",          dirk2,         2},
        {"qin_zhang",      qinZhang,      2},
        {"rk3",            rk3,           3},
        {"dirk3",          dirk3,         3},
        {"rk4",            rk4,           4}
    };

    const std::vector<Scalar> timeStepSizes = {0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625};

    Dumux::Json::JsonTree out;
    out["methods"] = Dumux::Json::JsonTree::object();
    for (const auto& info : allMethods)
    {
        auto& node = out["methods"][info.id];
        node["name"] = displayName(info.id, info.method);
        node["order"] = info.order;

        std::vector<std::array<Scalar, 2>> convergence;
        convergence.reserve(timeStepSizes.size());
        for (const auto dt : timeStepSizes)
            convergence.push_back({dt, integrateToLin(info.method, dt)});
        node["convergence"] = convergence;
    }

    const std::string outputFile = "test_timestepmethods_data.json";
    std::ofstream output(outputFile);
    output << std::setw(2) << out << std::endl;
    std::cout << "Wrote convergence data to " << outputFile << std::endl;

    return 0;
}
