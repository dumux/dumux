// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief A time stepper performing a single time step of a transient simulation
 */
#ifndef DUMUX_TIMESTEPPING_MULTISTAGE_TIMESTEPPER_HH
#define DUMUX_TIMESTEPPING_MULTISTAGE_TIMESTEPPER_HH

#include <memory>
#include <vector>
#include <cmath>
#include <iostream>

#include <dumux/io/format.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/variablesbackend.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/experimental/timestepping/timelevel.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>
#include <dumux/experimental/timestepping/multistagemethods.hh>

namespace Dumux::Experimental {

//! Data object for the parameters of a given stage
template<class Scalar>
class MultiStageParams
{
    struct Params {
        Scalar alpha, betaDt, timeAtStage, dtFraction;
        bool skipTemporal, skipSpatial;
    };
public:
    //! Extract params for stage i from method m
    MultiStageParams(const MultiStageMethod<Scalar>& m, std::size_t i, const Scalar t, const Scalar dt)
    : size_(i+1)
    {
        params_.resize(size_);
        for (std::size_t k = 0; k < size_; ++k)
        {
            auto& p = params_[k];
            p.alpha = m.temporalWeight(i, k);
            p.betaDt = m.spatialWeight(i, k)*dt;
            p.timeAtStage = t + m.timeStepWeight(k)*dt;
            p.dtFraction = m.timeStepWeight(k);

            // decide on the dimensionless tableau coefficients; betaDt carries the
            // time step and is small for an active stage when dt is small
            using std::abs;
            p.skipTemporal = (abs(m.temporalWeight(i, k)) < 1e-6);
            p.skipSpatial = (abs(m.spatialWeight(i, k)) < 1e-6);
        }
    }

    std::size_t size () const
    { return size_; }

    //! weights of the temporal operator residual (\f$ \alpha_{ik} \f$)
    Scalar temporalWeight (std::size_t k) const
    { return params_[k].alpha; }

    //! weights of the spatial operator residual (\f$ \beta_{ik} \f$)
    Scalar spatialWeight (std::size_t k) const
    { return params_[k].betaDt; }

    //! the time at which we have to evaluate the operators
    Scalar timeAtStage (std::size_t k) const
    { return params_[k].timeAtStage; }

    //! the fraction of a time step corresponding to the k-th stage
    Scalar timeStepFraction (std::size_t k) const
    { return params_[k].dtFraction; }

    //! If \f$ \alpha_{ik} = 0\f$
    bool skipTemporal (std::size_t k) const
    { return params_[k].skipTemporal; }

    //! If \f$ \beta_{ik} = 0\f$
    bool skipSpatial (std::size_t k) const
    { return params_[k].skipSpatial; }

private:
    std::size_t size_;
    std::vector<Params> params_;
};

/*!
 * \brief Time stepping with a multi-stage method
 * \note We limit ourselves to "diagonally" implicit multi-stage methods where solving
 *       a stage can only depend on the values of the same stage and stages before
 *       but not future stages (which would require solving larger linear systems)
 * \note The stepper evaluates the storage and spatial operators at the solution of each
 *       stage and hands their weighted sum to the assembler (see ImplicitSolveParams).
 *       The residual of each stage is scaled by \f$ 1/\Delta t \f$, i.e. its spatial
 *       part has the units of the spatial operator.
 */
template<class PDESolver, class Scalar = double>
class MultiStageTimeStepper
{
    using Variables = typename PDESolver::Variables;
    using Assembler = typename PDESolver::Assembler;
    using Residual = typename Assembler::ResidualType;
    using StageParams = MultiStageParams<Scalar>;
    using SolveParams = ImplicitSolveParams<Scalar, Residual>;
    using Backend = VariablesBackend<Variables>;

public:

    /*!
     * \brief The constructor
     * \param pdeSolver Solver class for solving a PDE in each stage
     * \param msMethod The multi-stage method which is to be used for time integration
     * \param paramGroup A parameter group in which we look for parameters
     */
    MultiStageTimeStepper(std::shared_ptr<PDESolver> pdeSolver,
                          std::shared_ptr<const MultiStageMethod<Scalar>> msMethod,
                          const std::string& paramGroup = "")
    : pdeSolver_(pdeSolver)
    , msMethod_(msMethod)
    {
        initParams_(paramGroup);

        if (pdeSolver_->verbosity() >= 1)
            std::cout << "Initialize time stepper with method " << msMethod_->name()
                    << Fmt::format(" ({} stage{})", msMethod_->numStages(), (msMethod_->numStages() > 1 ? "s" : ""))
                    << std::endl;
    }

    /*!
     * \brief Advance one time step of the given time loop
     * \param vars The variables object at the current time level.
     * \param t The current time level
     * \param dt The time step size to be used
     * \note We expect the time level in vars to correspond to the given time `t`
     */
    void step(Variables& vars, const Scalar t, const Scalar dt)
    {
        if (!step_(vars, t, dt))
            DUNE_THROW(NumericalProblem, "Solver did not converge!");
    }

    /*!
     * \brief Try to advance one time step
     * \param vars The variables object at the current time level
     * \param t The current time level
     * \param dt The time step size to be used
     * \return whether the solver converged in all stages
     * \note If the solver did not converge, vars may be in an intermediate state.
     */
    bool tryStep(Variables& vars, const Scalar t, const Scalar dt)
    { return step_(vars, t, dt); }

    /*!
     * \brief Advance one time step of the given time loop (adaptive time stepping on solver failure)
     * \param vars The variables object at the current time level.
     * \param timeLoop An instance of a time loop
     * \note We expect the time level in vars to correspond to the given time `t`
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

private:
    bool step_(Variables& vars, const Scalar t, const Scalar dt)
    {
        auto& assembler = pdeSolver_->assembler();
        const auto numStages = msMethod_->numStages();
        storage_.resize(numStages);
        spatial_.resize(numStages);

        assembler.assembleOperators(vars, TimeLevel<Scalar>{t, t, msMethod_->timeStepWeight(0)},
                                    {.storage = &storage_[0], .spatial = &spatial_[0]});

        for (std::size_t stageIdx = 1; stageIdx <= numStages; ++stageIdx)
        {
            const StageParams params(*msMethod_, stageIdx, t, dt);
            const TimeLevel<Scalar> timeLevel{params.timeAtStage(stageIdx), t, params.timeStepFraction(stageIdx)};

            auto constant = std::make_shared<Residual>(storage_[0]);
            *constant = 0.0;
            for (std::size_t k = 0; k < stageIdx; ++k)
            {
                if (!params.skipTemporal(k))
                    Detail::TimeStepping::axpy(*constant, msMethod_->temporalWeight(stageIdx, k)/dt, storage_[k]);
                if (!params.skipSpatial(k))
                    Detail::TimeStepping::axpy(*constant, msMethod_->spatialWeight(stageIdx, k), spatial_[k]);
            }

            OperatorWeights<Scalar> weights;
            if (!params.skipTemporal(stageIdx))
                weights.storage = msMethod_->temporalWeight(stageIdx, stageIdx)/dt;
            if (!params.skipSpatial(stageIdx))
                weights.spatial = msMethod_->spatialWeight(stageIdx, stageIdx);

            assembler.setImplicitSolve(vars, SolveParams{timeLevel, weights, std::move(constant)});

            if (!pdeSolver_->apply(vars))
                return false;

            if (stageIdx < numStages)
                assembler.assembleOperators(vars, timeLevel,
                                            {.storage = &storage_[stageIdx], .spatial = &spatial_[stageIdx]});
        }

        return true;
    }

    void initParams_(const std::string& group = "")
    {
        maxTimeStepDivisions_ = getParamFromGroup<std::size_t>(group, "TimeStepper.MaxTimeStepDivisions", 10);
        retryTimeStepReductionFactor_ = getParamFromGroup<Scalar>(group, "TimeStepper.RetryTimeStepReductionFactor", 0.5);
    }

    std::shared_ptr<PDESolver> pdeSolver_;
    std::shared_ptr<const MultiStageMethod<Scalar>> msMethod_;

    std::vector<Residual> storage_; //!< storage operator at the solutions of the stages
    std::vector<Residual> spatial_; //!< spatial operator at the solutions of the stages

    std::size_t maxTimeStepDivisions_;
    Scalar retryTimeStepReductionFactor_;
};

} // end namespace Dumux::Experimental

#endif
