// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Experimental
 * \brief The generalized-alpha method for problems with second time derivatives
 */
#ifndef DUMUX_TIMESTEPPING_GENERALIZED_ALPHA_HH
#define DUMUX_TIMESTEPPING_GENERALIZED_ALPHA_HH

#include <memory>
#include <string>
#include <utility>
#include <optional>
#include <iostream>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/io/format.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/variablesbackend.hh>
#include <dumux/common/timeloop.hh>

#include <dumux/experimental/timestepping/timelevel.hh>
#include <dumux/experimental/timestepping/implicitsolve.hh>
#include <dumux/experimental/timestepping/timeintegrationmethod.hh>

namespace Dumux::Experimental {

/*!
 * \ingroup Experimental
 * \brief Parameters of the generalized-alpha method of Chung and Hulbert (1993) \cite ChungHulbert1993
 *
 * For the semi-discrete problem
 *
 * \f[
 *   \frac{\partial^2 I(x, t)}{\partial t^2} + \frac{\partial M(x, t)}{\partial t} + R(x, t) = 0
 * \f]
 *
 * (see ImplicitSolveParams), the method advances the values of the inertia operator \f$ I \f$
 * and the storage operator \f$ M \f$ together with their first and second time derivatives
 * with the Newmark relations
 *
 * \f[
 *   X^{n+1} = X^n + \Delta t \dot X^n + \Delta t^2 \left[ \left(\tfrac{1}{2} - \beta\right) \ddot X^n + \beta \ddot X^{n+1} \right], \quad
 *   \dot X^{n+1} = \dot X^n + \Delta t \left[ (1 - \gamma) \ddot X^n + \gamma \ddot X^{n+1} \right]
 * \f]
 *
 * for \f$ X \in \{I, M\} \f$ and solves
 *
 * \f[
 *   (1 - \alpha_m) \ddot I^{n+1} + \alpha_m \ddot I^n + (1 - \alpha_f) \dot M^{n+1} + \alpha_f \dot M^n
 *     + (1 - \alpha_f) R(x^{n+1}, t^{n+1}) + \alpha_f R(x^n, t^n) = 0.
 * \f]
 *
 * For linear \f$ I \f$ and \f$ M \f$ this is the method of Chung and Hulbert for
 * \f$ \mathbf{M}\ddot{\mathbf{u}} + \mathbf{C}\dot{\mathbf{u}} + \mathbf{r}(\mathbf{u}) = 0 \f$,
 * with the spatial operator averaged between the time levels. With
 * \f$ \gamma = \tfrac{1}{2} - \alpha_m + \alpha_f \f$ the method is second order accurate, and with
 * \f$ \beta = \tfrac{1}{4}(1 - \alpha_m + \alpha_f)^2 \f$ and \f$ \alpha_m \leq \alpha_f \leq \tfrac12 \f$
 * unconditionally stable for linear problems. Equations without inertia and storage
 * (constraints) are damped by the factor \f$ -\alpha_f/(1-\alpha_f) \f$ per step.
 */
template<class Scalar>
class GeneralizedAlphaMethod : public TimeIntegrationMethod<Scalar>
{
public:
    GeneralizedAlphaMethod(const Scalar alphaM, const Scalar alphaF,
                           const Scalar beta, const Scalar gamma,
                           const std::string& name)
    : alphaM_(alphaM), alphaF_(alphaF), beta_(beta), gamma_(gamma), name_(name)
    {
        if (!(beta_ > 0.0))
            DUNE_THROW(Dune::InvalidStateException, "The generalized-alpha method requires beta > 0");
        if (!(alphaM_ < 1.0 && alphaF_ < 1.0))
            DUNE_THROW(Dune::InvalidStateException, "The generalized-alpha method requires alpha_m < 1 and alpha_f < 1");
    }

    bool implicit () const final
    { return true; }

    std::string name () const final
    { return name_; }

    //! the weight \f$ \alpha_m \f$ of the previous time level in the inertia term
    Scalar alphaM () const { return alphaM_; }

    //! the weight \f$ \alpha_f \f$ of the previous time level in the storage and spatial terms
    Scalar alphaF () const { return alphaF_; }

    //! the Newmark parameter \f$ \beta \f$
    Scalar beta () const { return beta_; }

    //! the Newmark parameter \f$ \gamma \f$
    Scalar gamma () const { return gamma_; }

protected:
    //! Second order accurate and unconditionally stable parameters for given \f$ \alpha_m, \alpha_f \f$
    static GeneralizedAlphaMethod fromAlphas_(const Scalar alphaM, const Scalar alphaF, const std::string& name)
    {
        const Scalar gamma = 0.5 - alphaM + alphaF;
        const Scalar beta = 0.25*(1.0 - alphaM + alphaF)*(1.0 - alphaM + alphaF);
        return GeneralizedAlphaMethod(alphaM, alphaF, beta, gamma, name);
    }

    static Scalar checkedSpectralRadius_(const Scalar rhoInf, const Scalar min)
    {
        if (!(rhoInf >= min && rhoInf <= 1.0))
            DUNE_THROW(Dune::InvalidStateException,
                "The spectral radius at infinity has to be in [" << min << ", 1], not " << rhoInf);
        return rhoInf;
    }

private:
    Scalar alphaM_, alphaF_, beta_, gamma_;
    std::string name_;
};

//! Parameter sets of the generalized-alpha method
namespace GeneralizedAlpha {

/*!
 * \brief The generalized-alpha method of Chung and Hulbert (1993) \cite ChungHulbert1993
 * \note With the spectral radius \f$ \rho_\infty \in [0, 1] \f$ for \f$ \omega\Delta t \to \infty \f$:
 *       \f$ \alpha_m = (2\rho_\infty - 1)/(\rho_\infty + 1) \f$, \f$ \alpha_f = \rho_\infty/(\rho_\infty + 1) \f$.
 *       The method has the least low-frequency dissipation for given high-frequency dissipation.
 */
template<class Scalar>
class ChungHulbert : public GeneralizedAlphaMethod<Scalar>
{
    using ParentType = GeneralizedAlphaMethod<Scalar>;
public:
    explicit ChungHulbert(const Scalar rhoInf)
    : ParentType(ParentType::fromAlphas_(
        (2.0*ParentType::checkedSpectralRadius_(rhoInf, 0.0) - 1.0)/(rhoInf + 1.0),
        rhoInf/(rhoInf + 1.0),
        Fmt::format("generalized-alpha (Chung-Hulbert, rho_inf = {})", rhoInf)))
    {}
};

/*!
 * \brief The method of Hilber, Hughes and Taylor (1977) \cite HilberHughesTaylor1977
 * \note With the spectral radius \f$ \rho_\infty \in [\tfrac12, 1] \f$:
 *       \f$ \alpha_m = 0 \f$, \f$ \alpha_f = (1 - \rho_\infty)/(1 + \rho_\infty) \f$.
 */
template<class Scalar>
class HilberHughesTaylor : public GeneralizedAlphaMethod<Scalar>
{
    using ParentType = GeneralizedAlphaMethod<Scalar>;
public:
    explicit HilberHughesTaylor(const Scalar rhoInf)
    : ParentType(ParentType::fromAlphas_(
        0.0,
        (1.0 - ParentType::checkedSpectralRadius_(rhoInf, 0.5))/(1.0 + rhoInf),
        Fmt::format("Hilber-Hughes-Taylor (rho_inf = {})", rhoInf)))
    {}
};

/*!
 * \brief The method of Wood, Bossak and Zienkiewicz (1980) \cite WoodBossakZienkiewicz1980
 * \note With the spectral radius \f$ \rho_\infty \in [0, 1] \f$:
 *       \f$ \alpha_m = (\rho_\infty - 1)/(\rho_\infty + 1) \f$, \f$ \alpha_f = 0 \f$.
 */
template<class Scalar>
class WoodBossakZienkiewicz : public GeneralizedAlphaMethod<Scalar>
{
    using ParentType = GeneralizedAlphaMethod<Scalar>;
public:
    explicit WoodBossakZienkiewicz(const Scalar rhoInf)
    : ParentType(ParentType::fromAlphas_(
        (ParentType::checkedSpectralRadius_(rhoInf, 0.0) - 1.0)/(rhoInf + 1.0),
        0.0,
        Fmt::format("Wood-Bossak-Zienkiewicz (rho_inf = {})", rhoInf)))
    {}
};

/*!
 * \brief The Newmark method (1959) \cite Newmark1959
 * \note \f$ \alpha_m = \alpha_f = 0 \f$. The default parameters \f$ \beta = \tfrac14 \f$,
 *       \f$ \gamma = \tfrac12 \f$ (average acceleration) give a second order, unconditionally
 *       stable method without numerical dissipation.
 */
template<class Scalar>
class Newmark : public GeneralizedAlphaMethod<Scalar>
{
    using ParentType = GeneralizedAlphaMethod<Scalar>;
public:
    explicit Newmark(const Scalar beta = 0.25, const Scalar gamma = 0.5)
    : ParentType(0.0, 0.0, beta, gamma, Fmt::format("Newmark (beta = {}, gamma = {})", beta, gamma))
    {}
};

} // end namespace GeneralizedAlpha

/*!
 * \ingroup Experimental
 * \brief Time stepping with the generalized-alpha method (see GeneralizedAlphaMethod)
 *
 * The stepper keeps the values of the inertia and storage operators with their first and
 * second time derivatives and the value of the spatial operator at the current time level.
 * The initial time derivatives are computed from given time derivatives of the solution with
 * the Jacobians of the operators (exact for operators linear in the solution). The residual of
 * each implicit solve has the units of the spatial operator.
 *
 * Problems without inertia can be integrated as well: the storage operator then follows the
 * Newmark relations. With the Chung-Hulbert parameters, the method then has the eigenvalues of the
 * first-order generalized-alpha method of Jansen, Whiting and Hulbert (2000) \cite JansenWhitingHulbert2000
 * with the same \f$ \rho_\infty \f$ and the additional eigenvalue \f$ -\alpha_f/(1-\alpha_f) \f$ of the
 * residual \f$ \dot M + R \f$.
 */
template<class PDESolver, class Scalar = double>
class GeneralizedAlphaTimeStepper
{
    using Variables = typename PDESolver::Variables;
    using Assembler = typename PDESolver::Assembler;
    using Residual = typename Assembler::ResidualType;
    using SolveParams = ImplicitSolveParams<Scalar, Residual>;
    using Backend = VariablesBackend<Variables>;
    using DofVector = std::decay_t<decltype(Backend::dofs(std::declval<const Variables&>()))>;

    static constexpr bool withInertia = []{
        if constexpr (requires { Assembler::hasInertia(); })
            return Assembler::hasInertia();
        else
            return false;
    }();

    //! value, first and second time derivative of an operator
    struct OperatorHistory
    {
        Residual value, rate, secondRate;
    };

public:
    /*!
     * \brief The constructor
     * \param pdeSolver Solver class for solving a PDE in each time step
     * \param method The parameters of the generalized-alpha method
     * \param paramGroup A parameter group in which we look for parameters
     */
    GeneralizedAlphaTimeStepper(std::shared_ptr<PDESolver> pdeSolver,
                                std::shared_ptr<const GeneralizedAlphaMethod<Scalar>> method,
                                const std::string& paramGroup = "")
    : pdeSolver_(pdeSolver)
    , method_(method)
    {
        maxTimeStepDivisions_ = getParamFromGroup<std::size_t>(paramGroup, "TimeStepper.MaxTimeStepDivisions", 10);
        retryTimeStepReductionFactor_ = getParamFromGroup<Scalar>(paramGroup, "TimeStepper.RetryTimeStepReductionFactor", 0.5);

        if (pdeSolver_->verbosity() >= 1)
            std::cout << "Initialize time stepper with method " << method_->name() << std::endl;
    }

    /*!
     * \brief Start the time integration at time t0
     * \param vars The variables object at t0
     * \param t0 The initial time
     * \note With inertia, the integration starts at rest (vanishing time derivatives of the operators).
     *       Without inertia, the storage operator starts with the rate \f$ -R(x_0, t_0) \f$ and a
     *       vanishing second time derivative. This rate is consistent if the initial solution
     *       satisfies the equations without storage term.
     */
    void initialize(Variables& vars, const Scalar t0)
    { initialize_(vars, t0, nullptr, nullptr); }

    /*!
     * \brief Start the time integration at time t0 with given time derivatives of the solution
     * \param vars The variables object at t0
     * \param t0 The initial time
     * \param velocity The first time derivative of the solution at t0
     * \param acceleration The second time derivative of the solution at t0
     * \note The time derivatives of the operators are their Jacobians times the time derivatives
     *       of the solution, which neglects second derivatives of nonlinear operators.
     */
    void initialize(Variables& vars, const Scalar t0, const DofVector& velocity, const DofVector& acceleration)
    { initialize_(vars, t0, &velocity, &acceleration); }

    /*!
     * \brief Advance one time step
     * \param vars The variables object at the current time level
     * \param t The current time level
     * \param dt The time step size to be used
     * \note Without a preceding call of initialize(), the integration starts at t as with initialize(vars, t).
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
     * \return whether the solver converged
     * \note If the solver did not converge, vars may be in an intermediate state.
     */
    bool tryStep(Variables& vars, const Scalar t, const Scalar dt)
    { return step_(vars, t, dt); }

    //! The value, first and second time derivative of the inertia operator at the current time level
    const auto& inertiaHistory () const
    { return inertia_; }

    //! The value, first and second time derivative of the storage operator at the current time level
    const auto& storageHistory () const
    { return storage_; }

private:
    void initialize_(Variables& vars, const Scalar t0, const DofVector* velocity, const DofVector* acceleration)
    {
        auto& assembler = pdeSolver_->assembler();
        assembler.assembleOperators(vars, TimeLevel<Scalar>{t0},
                                    {.inertia = withInertia ? &inertia_.value : nullptr,
                                     .storage = &storage_.value, .spatial = &spatial_});

        for (auto* history : {&inertia_, &storage_})
        {
            history->rate = storage_.value; history->rate = 0.0;
            history->secondRate = storage_.value; history->secondRate = 0.0;
        }
        if constexpr (!withInertia)
        {
            inertia_.value = storage_.value;
            inertia_.value = 0.0;
        }

        if (velocity && acceleration)
        {
            const auto applyJacobian = [&](const OperatorWeights<Scalar>& weights, OperatorHistory& history)
            {
                assembler.setImplicitSolve(vars, SolveParams{TimeLevel<Scalar>{t0}, weights, nullptr});
                assembler.assembleJacobianAndResidual(vars);
                const auto& jacobian = assembler.jacobian();
                if constexpr (requires { jacobian.mv(*velocity, history.rate); })
                {
                    jacobian.mv(*velocity, history.rate);
                    jacobian.mv(*acceleration, history.secondRate);
                }
                else
                {
                    history.rate = jacobian*(*velocity);
                    history.secondRate = jacobian*(*acceleration);
                }
            };

            if constexpr (withInertia)
                applyJacobian(OperatorWeights<Scalar>{1.0, 0.0, 0.0}, inertia_);
            applyJacobian(OperatorWeights<Scalar>{0.0, 1.0, 0.0}, storage_);
        }
        else if constexpr (!withInertia)
        {
            // the equation dM/dt + R = 0 at t0; the second time derivative does not affect the order
            storage_.rate = spatial_;
            storage_.rate *= -1.0;
        }

        time_ = t0;
    }

    bool step_(Variables& vars, const Scalar t, const Scalar dt)
    {
        if (!time_)
            initialize(vars, t);
        else if (!Dune::FloatCmp::eq(*time_, t))
            DUNE_THROW(Dune::InvalidStateException,
                "The generalized-alpha time stepper is at time " << *time_ << " and cannot step from time " << t
                << ". Call initialize() to restart the time integration.");

        const auto am = method_->alphaM();
        const auto af = method_->alphaF();
        const auto beta = method_->beta();
        const auto gamma = method_->gamma();
        using Detail::TimeStepping::axpy;

        // the second time derivative at the new time level is X(x)/(beta dt^2) + rest
        const auto secondRateRest = [&](const OperatorHistory& h)
        {
            Residual rest = h.value;
            rest *= -1.0/(beta*dt*dt);
            axpy(rest, -1.0/(beta*dt), h.rate);
            axpy(rest, -(0.5 - beta)/beta, h.secondRate);
            return rest;
        };

        auto constant = std::make_shared<Residual>(storage_.value);
        *constant = 0.0;

        if constexpr (withInertia)
        {
            axpy(*constant, 1.0 - am, secondRateRest(inertia_));
            axpy(*constant, am, inertia_.secondRate);
        }

        // first time derivative of the storage at the new time level is gamma dt/(beta dt^2) M(x) + rest
        auto storageRateRest = storage_.rate;
        axpy(storageRateRest, dt*(1.0 - gamma), storage_.secondRate);
        axpy(storageRateRest, gamma*dt, secondRateRest(storage_));
        axpy(*constant, 1.0 - af, storageRateRest);
        axpy(*constant, af, storage_.rate);

        if (af != 0.0)
            axpy(*constant, af, spatial_);

        OperatorWeights<Scalar> weights;
        weights.inertia = withInertia ? (1.0 - am)/(beta*dt*dt) : 0.0;
        weights.storage = (1.0 - af)*gamma/(beta*dt);
        weights.spatial = 1.0 - af;

        auto& assembler = pdeSolver_->assembler();
        assembler.setImplicitSolve(vars, SolveParams{TimeLevel<Scalar>{t + dt, t, 1.0}, weights, std::move(constant)});
        if (!pdeSolver_->apply(vars))
            return false;

        // advance the operator histories to the new time level
        OperatorHistory newInertia, newStorage;
        Residual newSpatial;
        assembler.assembleOperators(vars, TimeLevel<Scalar>{t + dt, t, 1.0},
                                    {.inertia = withInertia ? &newInertia.value : nullptr,
                                     .storage = &newStorage.value, .spatial = &newSpatial});

        const auto advance = [&](OperatorHistory& h, OperatorHistory& next)
        {
            next.secondRate = secondRateRest(h);
            axpy(next.secondRate, 1.0/(beta*dt*dt), next.value);
            next.rate = h.rate;
            axpy(next.rate, dt*(1.0 - gamma), h.secondRate);
            axpy(next.rate, dt*gamma, next.secondRate);
            h = std::move(next);
        };

        if constexpr (withInertia)
            advance(inertia_, newInertia);
        advance(storage_, newStorage);
        spatial_ = std::move(newSpatial);

        time_ = t + dt;
        return true;
    }

    std::shared_ptr<PDESolver> pdeSolver_;
    std::shared_ptr<const GeneralizedAlphaMethod<Scalar>> method_;

    OperatorHistory inertia_, storage_;
    Residual spatial_;
    std::optional<Scalar> time_;

    std::size_t maxTimeStepDivisions_;
    Scalar retryTimeStepReductionFactor_;
};

} // end namespace Dumux::Experimental

#endif
