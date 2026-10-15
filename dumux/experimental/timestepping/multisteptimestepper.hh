// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief A time stepper for linear multi-step methods
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTEP_TIMESTEPPER_HH
#define DUMUX_TIMESTEPPING_MULTISTEP_TIMESTEPPER_HH

#include <deque>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <vector>
#include <iostream>

#include <dune/common/float_cmp.hh>

#include <dumux/io/format.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/variablesbackend.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/experimental/timestepping/timelevel.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>
#include <dumux/experimental/timestepping/multistepmethods.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>
#include <dumux/experimental/timestepping/multistagetimestepper.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \brief Time stepping with a linear multi-step method (see MultiStepMethod)
 *
 * The stepper evaluates the storage operator (and the spatial operator if the method
 * requires it) at the beginning of each time step and keeps the evaluations of as many
 * previous time levels as the method uses. The residual of each implicit solve is scaled
 * by \f$ 1/\Delta t \f$, i.e. its spatial part has the units of the spatial operator.
 *
 * While fewer previous time levels are available than the method uses, i.e. at the
 * beginning of a simulation and after resetHistory(), the steps are taken with a startup
 * method if one is given, and otherwise with the number of levels the method selects for
 * the available history. For BDF methods of order three and higher, the default startup
 * method is the fourth order DIRK method of Hairer and Wanner, which keeps the global
 * order of BDF3 to BDF5. Starting BDF1 and BDF2 with lower BDF orders keeps their order.
 *
 * The history is discarded automatically if a step does not start at the end of the
 * previous step, e.g. after restarting at a different time.
 */
template<class PDESolver, class Scalar = double>
class MultiStepTimeStepper
{
    using Variables = typename PDESolver::Variables;
    using Assembler = typename PDESolver::Assembler;
    using Residual = typename Assembler::ResidualType;
    using SolveParams = ImplicitSolveParams<Scalar, Residual>;
    using Backend = VariablesBackend<Variables>;

    struct TimeLevelData
    {
        Scalar time;
        Residual storage;
        Residual spatial;
    };

public:
    /*!
     * \brief The constructor with the default startup method
     * \param pdeSolver Solver class for solving a PDE in each time step
     * \param method The multi-step method which is to be used for time integration
     * \param paramGroup A parameter group in which we look for parameters
     */
    MultiStepTimeStepper(std::shared_ptr<PDESolver> pdeSolver,
                         std::shared_ptr<const MultiStepMethod<Scalar>> method,
                         const std::string& paramGroup = "")
    : MultiStepTimeStepper(pdeSolver, method, defaultStartupMethod_(*method), paramGroup)
    {}

    /*!
     * \brief The constructor
     * \param pdeSolver Solver class for solving a PDE in each time step
     * \param method The multi-step method which is to be used for time integration
     * \param startupMethod The method for steps with fewer previous time levels than the method uses
     *        (nullptr to start with fewer time levels of the multi-step method)
     * \param paramGroup A parameter group in which we look for parameters
     */
    MultiStepTimeStepper(std::shared_ptr<PDESolver> pdeSolver,
                         std::shared_ptr<const MultiStepMethod<Scalar>> method,
                         std::shared_ptr<const MultiStageMethod<Scalar>> startupMethod,
                         const std::string& paramGroup = "")
    : pdeSolver_(pdeSolver)
    , method_(method)
    {
        static_assert(!assemblerHasInertia_(),
            "The multi-step time stepper integrates first time derivatives only, but the problem enables inertia");

        initParams_(paramGroup);

        if (startupMethod)
            startupStepper_ = std::make_unique<MultiStageTimeStepper<PDESolver, Scalar>>(pdeSolver, startupMethod, paramGroup);

        if (pdeSolver_->verbosity() >= 1)
            std::cout << "Initialize time stepper with method " << method_->name()
                      << (startupMethod ? " (startup with " + startupMethod->name() + ")" : std::string{})
                      << std::endl;
    }

    /*!
     * \brief Advance one time step
     * \param vars The variables object at the current time level
     * \param t The current time level
     * \param dt The time step size to be used
     */
    void step(Variables& vars, const Scalar t, const Scalar dt)
    {
        if (!step_(vars, t, dt))
            DUNE_THROW(NumericalProblem, "Solver did not converge!");
    }

    /*!
     * \brief Advance one time step of the given time loop (adaptive time stepping on solver failure)
     * \param vars The variables object at the current time level
     * \param timeLoop An instance of a time loop
     */
    void step(Variables& vars, TimeLoopBase<Scalar>& timeLoop)
    {
        const auto initialDofs = Backend::dofs(vars);

        for (std::size_t i = 0; i <= maxTimeStepDivisions_; ++i)
        {
            if (step_(vars, timeLoop.time(), timeLoop.timeStepSize()))
                return;

            if (i < maxTimeStepDivisions_)
            {
                Backend::update(vars, initialDofs);
                pdeSolver_->assembler().resetTimeStep(Backend::dofs(vars));

                if (pdeSolver_->verbosity() >= 1)
                {
                    const auto dt = timeLoop.timeStepSize();
                    std::cout << Fmt::format("Solver did not converge with dt = {} seconds. ", dt)
                              << Fmt::format("Retrying with time step of dt = {} seconds.\n", dt*retryTimeStepReductionFactor_);
                }

                timeLoop.setTimeStepSize(timeLoop.timeStepSize() * retryTimeStepReductionFactor_);
            }
        }

        DUNE_THROW(NumericalProblem,
            Fmt::format("Solver didn't converge after {} time-step divisions; dt = {}.\n",
                        maxTimeStepDivisions_, timeLoop.timeStepSize()));
    }

    /*!
     * \brief Try to advance one time step
     * \param vars The variables object at the current time level
     * \param t The current time level
     * \param dt The time step size to be used
     * \return whether the solver converged
     * \note If the solver did not converge, vars may be in an intermediate state. A repeated
     *       call from the same time level reuses the evaluations at that level, so vars has to
     *       be reset to the solution at that level before.
     */
    bool tryStep(Variables& vars, const Scalar t, const Scalar dt)
    { return step_(vars, t, dt); }

    /*!
     * \brief Discard the evaluations at previous time levels
     * \note Required whenever the solution changes between steps in a way the stepper does not see,
     *       e.g. after grid adaptation.
     */
    void resetHistory()
    {
        history_.clear();
        stepSizes_.clear();
        lastStepEnd_.reset();
        pending_ = false;
    }

    //! The number of previous time levels used in the last step (zero for a step with the startup method)
    std::size_t lastNumSteps() const
    { return lastNumSteps_; }

private:
    bool step_(Variables& vars, const Scalar t, const Scalar dt)
    {
        auto& assembler = pdeSolver_->assembler();

        const bool retry = pending_ && Dune::FloatCmp::eq(history_.front().time, t);
        if (!retry)
        {
            if (!lastStepEnd_ || !Dune::FloatCmp::eq(*lastStepEnd_, t))
                resetHistory();

            history_.emplace_front();
            auto& current = history_.front();
            current.time = t;
            assembler.assembleOperators(vars, TimeLevel<Scalar>{t, t, 0.0},
                                        {.storage = &current.storage,
                                         .spatial = method_->spatialHistory() ? &current.spatial : nullptr});
            pending_ = true;
        }

        bool converged = false;
        if (startupStepper_ && history_.size() < method_->numSteps())
        {
            lastNumSteps_ = 0;
            converged = startupStepper_->tryStep(vars, t, dt);
        }
        else
        {
            std::vector<Scalar> stepSizes{dt};
            stepSizes.insert(stepSizes.end(), stepSizes_.begin(), stepSizes_.end());
            lastNumSteps_ = method_->numStepsToUse(stepSizes);
            converged = solve_(vars, t, std::span<const Scalar>(stepSizes).first(lastNumSteps_));
        }

        if (converged)
        {
            pending_ = false;
            lastStepEnd_ = t + dt;
            stepSizes_.push_front(dt);
            while (history_.size() > method_->numSteps() - 1)
                history_.pop_back();
            while (stepSizes_.size() > method_->numSteps() - 1)
                stepSizes_.pop_back();
        }

        return converged;
    }

    //! Solve for the next time level with the given step sizes (one per used previous time level)
    bool solve_(Variables& vars, const Scalar t, std::span<const Scalar> stepSizes)
    {
        const auto k = stepSizes.size();
        const auto dt = stepSizes[0];
        std::vector<Scalar> temporal(k+1), spatial(k+1);
        method_->weights(stepSizes, temporal, spatial);

        auto constant = std::make_shared<Residual>(history_.front().storage);
        *constant = 0.0;
        for (std::size_t j = 1; j <= k; ++j)
        {
            Detail::TimeStepping::axpy(*constant, temporal[j]/dt, history_[j-1].storage);
            if (spatial[j] != 0.0)
                Detail::TimeStepping::axpy(*constant, spatial[j], history_[j-1].spatial);
        }

        OperatorWeights<Scalar> weights;
        weights.storage = temporal[0]/dt;
        weights.spatial = spatial[0];

        pdeSolver_->assembler().setImplicitSolve(vars, SolveParams{TimeLevel<Scalar>{t + dt, t, 1.0}, weights, std::move(constant)});
        return pdeSolver_->apply(vars);
    }

    static std::shared_ptr<const MultiStageMethod<Scalar>> defaultStartupMethod_(const MultiStepMethod<Scalar>& method)
    {
        if (method.numSteps() >= 3)
            return std::make_shared<MultiStage::DIRKFourthOrderHairerWanner<Scalar>>();
        return nullptr;
    }

    static constexpr bool assemblerHasInertia_()
    {
        if constexpr (requires { Assembler::hasInertia(); })
            return Assembler::hasInertia();
        else
            return false;
    }

    void initParams_(const std::string& group = "")
    {
        maxTimeStepDivisions_ = getParamFromGroup<std::size_t>(group, "TimeStepper.MaxTimeStepDivisions", 10);
        retryTimeStepReductionFactor_ = getParamFromGroup<Scalar>(group, "TimeStepper.RetryTimeStepReductionFactor", 0.5);
    }

    std::shared_ptr<PDESolver> pdeSolver_;
    std::shared_ptr<const MultiStepMethod<Scalar>> method_;
    std::unique_ptr<MultiStageTimeStepper<PDESolver, Scalar>> startupStepper_;

    std::deque<TimeLevelData> history_; //!< evaluations at the current and previous time levels, most recent first
    std::deque<Scalar> stepSizes_; //!< sizes of the previous steps, most recent first
    std::optional<Scalar> lastStepEnd_;
    bool pending_ = false; //!< whether the evaluations at the current time level belong to an unfinished step
    std::size_t lastNumSteps_ = 0;

    std::size_t maxTimeStepDivisions_;
    Scalar retryTimeStepReductionFactor_;
};

} // end namespace Dumux::Experimental

#endif
