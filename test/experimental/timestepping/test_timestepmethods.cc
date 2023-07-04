//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <utility>
#include <memory>
#include <array>

#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>

#include <dumux/io/format.hh>
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
    using StageParams = Experimental::MultiStageParams<Scalar>;

    void setLinearSystem() {}
    JacobianMatrix& jacobian() { return jac_; }
    ResidualType& residual() { return res_; }

    void assembleResidual(const Variables& vars)
    {
        if (stageParams_->size() != spatialResiduals_.size())
            DUNE_THROW(Dune::InvalidStateException, "Wrong number of residuals");

        // assemble residual components depending on current solution
        assembleResiduals_(vars,
            temporalResiduals_.back(),
            spatialResiduals_.back()
        );

        // assemble residual for current stage solver
        res_ = 0.0;
        for (std::size_t k = 0; k < stageParams_->size(); ++k)
        {
            if (!stageParams_->skipTemporal(k))
                res_ += stageParams_->temporalWeight(k)*temporalResiduals_[k];
            if (!stageParams_->skipSpatial(k))
                res_ += stageParams_->spatialWeight(k)*spatialResiduals_[k];
        }
    }

    void assembleJacobianAndResidual(const Variables& vars)
    {
        assembleResidual(vars);
        jac_ = 1.0;
    }

    void prepareStage(Variables& variables,
                      std::shared_ptr<const StageParams> params)
    {
        const auto curStage = params->size() - 1;
        std::cout << "Assembler is preparing stage " << curStage << std::endl;
        const auto prevStage = stageParams_ ? stageParams_->size() - 1 : 0;
        if (curStage != prevStage+1)
            DUNE_THROW(Dune::InvalidStateException,
                "Can only prepare stages consecutively (current stage: " << curStage
                << ", previous stage: " << prevStage << ")");

        stageParams_ = params;

        // in the first stage, also assemble the old residual
        if (curStage == 1)
        {
            // we update the time level of the given grid variables
            const auto t = params->timeAtStage(0);
            const auto prevT = params->timeAtStage(0);
            const auto dtFraction = params->timeStepFraction(0);
            variables.updateTime(Experimental::TimeLevel{t, prevT, dtFraction});

            // allocate memory for residuals of stage 0
            spatialResiduals_.emplace_back(0.0);
            temporalResiduals_.emplace_back(0.0);

            // assemble stage 0 residuals
            assembleResiduals_(variables,
                temporalResiduals_.back(),
                spatialResiduals_.back()
            );
        }

        // we update the time level of the given grid variables
        const auto t = params->timeAtStage(curStage);
        const auto prevT = params->timeAtStage(0);
        const auto dtFraction = params->timeStepFraction(curStage);
        variables.updateTime(Experimental::TimeLevel{t, prevT, dtFraction});

        // allocate memory for residuals of the upcoming stage
        spatialResiduals_.emplace_back(0.0);
        temporalResiduals_.emplace_back(0.0);
    }

    void clearStages()
    {
        std::cout << "Assembler is clearing stage data " << std::endl;
        spatialResiduals_.clear();
        temporalResiduals_.clear();
        stageParams_.reset();
    }

private:
    void assembleResiduals_(const Variables& variables,
                            ResidualType& temporal,
                            ResidualType& spatial) const
    {
        using std::exp;
        temporal = variables.dofs();
        spatial = -exp(variables.timeLevel().current());
    }

    ResidualType res_;
    JacobianMatrix jac_;
    std::vector<ResidualType> spatialResiduals_;
    std::vector<ResidualType> temporalResiduals_;
    std::shared_ptr<const StageParams> stageParams_;
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

    const auto testIntegration = [&] (auto method, Scalar expectedRelError)
    {
        std::cout << "\n-- Integration with " << method->name() << ":\n\n";
        SolutionVector x = 0.0;
        Variables vars(x);

        using TimeStepper = Experimental::MultiStageTimeStepper<NewtonSolver>;
        TimeStepper timeStepper(newtonSolver, method);

        const Scalar dt = 0.01;
        timeStepper.step(vars, /*time*/0.0, dt);

        const auto [abs, rel] = computeError(vars);
        std::cout << "\n"
                  << "-- Summary\n"
                  << "-- ===========================\n"
                  << "-- Exact solution:    " << Fmt::format("{:.4e}\n", exact(dt))
                  << "-- Discrete solution: " << Fmt::format("{:.4e}\n", vars.dofs())
                  << "-- Absolute error:    " << Fmt::format("{:.4e}\n", abs)
                  << "-- Relative error:    " << Fmt::format("{:.4e}", rel) << std::endl;

        if (Dune::FloatCmp::ne(rel, expectedRelError, 0.01))
            DUNE_THROW(Dune::InvalidStateException,
                       "Detected deviation from reference > 1% for " << method->name());
    };

    using namespace Experimental::MultiStage;
    testIntegration(std::make_shared<ExplicitEuler<Scalar>>(), 4.9917e-03);
    testIntegration(std::make_shared<ImplicitEuler<Scalar>>(), 5.0083e-03);
    testIntegration(std::make_shared<Theta<Scalar>>(0.5), 8.3333e-06);
    testIntegration(std::make_shared<DIRKThirdOrderAlexander<Scalar>>(), 7.9214e-09);
    testIntegration(std::make_shared<RungeKuttaExplicitFourthOrder<Scalar>>(), 3.4829e-12);

    return 0;
}
