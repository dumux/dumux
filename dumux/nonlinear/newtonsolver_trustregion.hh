// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Newton solver with a trust region approach
 */
#ifndef DUMUX_NEWTON_SOLVER_TRUST_REGION_HH
#define DUMUX_NEWTON_SOLVER_TRUST_REGION_HH

#include <string>
#include <dune/common/exceptions.hh>
#include <dumux/io/format.hh>
#include <dumux/common/parameters.hh>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \brief An implementation of a trust region Newton solver
 * \tparam Assembler the assembler
 * \tparam LinearSolver the linear solver
 * \tparam Comm the communication object used to communicate with all processes
 * \note The mathematical background is given in Chapter 4 of the book \cite NocedalWright2006
 */
template <class Assembler, class LinearSolver,
          class Reassembler = DefaultPartialReassembler,
          class Comm = Dune::Communication<Dune::MPIHelper::MPICommunicator> >
class NewtonSolverTrustRegion
: public NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>
{
    using ParentType = NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>;
    enum class FallbackStrategy { Newton, Dogleg, CauchyPoint };
    using Scalar = typename Assembler::Scalar;

protected:
    using Backend = typename ParentType::Backend;
    using SolutionVector = typename ParentType::SolutionVector;
    using ResidualVector = typename ParentType::ResidualVector;
    using LinearAlgebraNativeBackend = typename ParentType::LinearAlgebraNativeBackend;

public:
    using typename ParentType::Variables;
    using Communication = Comm;

    NewtonSolverTrustRegion(
        std::shared_ptr<Assembler> assembler,
        std::shared_ptr<LinearSolver> linearSolver,
        const Communication& comm = Dune::MPIHelper::getCommunication(),
        const std::string& paramGroup = "",
        const std::string& paramGroupName = "NewtonTrustRegion",
        int verbosity = 2
    )
    : ParentType(assembler, linearSolver, comm, paramGroup, paramGroupName, verbosity)
    {
        // Trust region radius parameters
        delta0_ = getParam<Scalar>("NewtonTrustRegion.TrustRegionDelta0", 1.0);
        deltaMin_ = getParam<Scalar>("NewtonTrustRegion.TrustRegionDeltaMin", 0.03);
        deltaMax_ = getParam<Scalar>("NewtonTrustRegion.TrustRegionDeltaMax", 3.0);
        incrementFactor_ = getParam<Scalar>("NewtonTrustRegion.TrustRegionIncrementFactor", 2.0);
        decrementFactor_ = getParam<Scalar>("NewtonTrustRegion.TrustRegionDecrementFactor", 0.25);
        incrementThreshold_ = getParam<Scalar>("NewtonTrustRegion.TrustRegionIncrementThreshold", 0.75);
        decrementThreshold_ = getParam<Scalar>("NewtonTrustRegion.TrustRegionDecrementThreshold", 0.25);
        eta_ = getParam<Scalar>("NewtonTrustRegion.TrustRegionAcceptanceThreshold", 0.25); // proposal acceptance threshold

        enableResidualCriterion_ = getParam<bool>("NewtonTrustRegion.EnableResidualCriterion", false);

        // Fallback strategy if the update is not within the trust region
        const auto strategy = getParam<std::string>("NewtonTrustRegion.TrustRegionFallbackStrategy", "Newton");
        if (strategy == "Newton")
            strategy_ = FallbackStrategy::Newton;
        else if (strategy == "Dogleg")
            strategy_ = FallbackStrategy::Dogleg;
        else if (strategy == "CauchyPoint")
            strategy_ = FallbackStrategy::CauchyPoint;
        else
            DUNE_THROW(Dune::Exception, "Unknown trust region fallback strategy: " << strategy
                        << ". Possible values are Newton, Dogleg, or CauchyPoint.");
    }

    bool enableResidualCriterion() const
    { return enableResidualCriterion_; }

    // The Newton implementation has a few hooks that we can override
    // The algorithm has the following structure where each step can be overridden:
    //
    //   newtonBegin
    //   while (newtonProceed(vars, converged))
    //       newtonBeginStep
    //       assembleLinearSystem
    //       solveLinearSystem
    //       newtonUpdate
    //       newtonEndStep
    //       newtonConverged
    //   newtonEnd
    //   if (!converged) newtonFail
    //   else newtonSucceed
    //

    void newtonBegin(Variables& vars) override
    {
        ParentType::newtonBegin(vars);

        // initialize temporary vectors
        g_ = LinearAlgebraNativeBackend::zeros(Backend::size(Backend::dofs(vars)));
        w_ = LinearAlgebraNativeBackend::zeros(Backend::size(Backend::dofs(vars)));
        deltaUU_ = LinearAlgebraNativeBackend::zeros(Backend::size(Backend::dofs(vars)));
    }

    void newtonBeginStep(const Variables& vars) override
    {
        ParentType::newtonBeginStep(vars);

        if (this->verbosity() >= 2)
        {
            const auto width = Fmt::formatted_size("{}", this->maxSteps_);
            std::cout << Fmt::format("\nNewtonTrustRegion begin iteration {:{}}\n", this->numSteps_+1, width);
        }
    }

    void assembleLinearSystem(const Variables& vars) override
    {
        // Assemble the residual and Jacobian if this is the first iteration
        // or if we need to reassemble due to a trust region update
        if (!alreadyAssembled_)
        {
            if (this->verbosity() >= 2)
                std::cout << "-- assemble residual and Jacobian" << std::endl;
            ParentType::assembleLinearSystem(vars); // computes the residual and Jacobian
        }
        else if (this->verbosity() >= 2)
            std::cout << "-- skipping assembly of residual and Jacobian" << std::endl;
    }

    // Main hook to update the trust region radius and decide if we accept the proposal step or not
    void newtonUpdate(Variables& vars, const SolutionVector& uLastIter, const ResidualVector& deltaU) override
    {
        // compute the norm of the proposed Newton step
        const auto newtonStepNorm = this->linearSolver().norm(deltaU);

        // the proposed step goes outside the trust region, invoke the fallback strategy
        if (newtonStepNorm > delta_ || newtonStepNorm == 0.0)
        {
            if (this->verbosity() >= 2)
                std::cout << Fmt::format("-- update outside trust region: ||Δu|| > Δ ({:.2e} > {:.2e})\n", newtonStepNorm, delta_);

            const auto fallback = newtonStepNorm > 0.0 ? strategy_ : FallbackStrategy::CauchyPoint;
            switch (fallback)
            {
                // take a Newton step with step size equal to the trust region radius
                case FallbackStrategy::Newton:
                    if (this->verbosity() >= 3)
                        std::cout << "-- fallback to restricted Newton update" << std::endl;
                    deltaU_ = deltaU;
                    alphaUpdate_ = delta_/newtonStepNorm;
                    break;

                // compute and use the Cauchy point as update (scaled gradient descent step)
                case FallbackStrategy::CauchyPoint:
                    if (this->verbosity() >= 3)
                        std::cout << "-- fallback to Cauchy point update" << std::endl;
                    computeCauchyPointUpdate_();
                    break;

                // take a combination of the Newton and Cauchy point steps: the Dogleg step
                case FallbackStrategy::Dogleg:
                    if (this->verbosity() >= 3)
                        std::cout << "-- fallback to Dogleg update" << std::endl;

                    // compute g^T B g where B = J^T J and g = J^T r
                    maybeAssembleGradientEtc_();

                    // Take the Cauchy point update
                    if (gTBg_ <= 0.0)
                        computeCauchyPointUpdate_();

                    // Dogleg step differentiates two cases (we already know that Newton step is outside the trust region):
                    // (A) Cauchy point is outside the trust region -> take restricted Cauchy step
                    // (B) Cauchy point is inside the trust region -> evaluate second leg of the dogleg step
                    else
                    {
                        // Cauchy point is outside the trust region
                        const auto alphaDogleg = gNorm_*gNorm_/gTBg_;
                        if (gNorm_*alphaDogleg >= delta_)
                        {
                            // take reduced gradient descent step predicted by the model
                            alphaUpdate_ = 1.0*delta_/gNorm_;
                            deltaU_ = g_;

                            if (this->verbosity() >= 3)
                                std::cout << Fmt::format("  -- Cauchy point outside trust region take restricted Cauchy step ||Δu|| = {:.2e}", alphaUpdate_*gNorm_);
                        }

                        // Cauchy point is inside the trust region: evaluate intersection of second leg
                        else
                        {
                            // Solve ||p_U + λ(p_B - p_U)||^2 - Δ^2 = 0 (see \cite NocedalWright2006 Ch.4/p.75)
                            // for solution λ = (τ - 1), τ in [1, 2] (τ for second leg)
                            // by bringing it in form a*λ^2 + b*λ + c = 0, where:
                            // a = ||p_B - p_U||^2, b = 2*p_U^T(p_B - p_U), c = ||p_U||^2 - Δ^2.
                            computeCauchyPointUpdate_();
                            const auto cauchyStepNorm = alphaUpdate_*gNorm_;
                            deltaUU_ = deltaU; // full Newton step (-p_B)
                            deltaU_ = 0.0; // scaled Cauchy point update (-p_U)
                            Backend::axpy(alphaDogleg, g_, deltaU_); // deltaU_ = ||g||^2/||g^T J^T J g|| * g
                            deltaUU_ -= deltaU_; // -(p_B - p_U)
                            const auto normDiff = this->linearSolver().norm(deltaUU_);  // ||(p_B - p_U)||

                            const auto a = normDiff*normDiff; // ||p_B - p_U||^2
                            const auto b = 2*(deltaUU_.dot(deltaU_)); // 2*p_U^T(p_B - p_U)
                            const auto c = cauchyStepNorm*cauchyStepNorm - delta_*delta_; // ||p_U||^2 - Δ^2
                            const auto roots = quadraticRoots_(a, b, c);

                            const auto tau = roots[0] >= 0.0 && roots[0] <= 1.0 ? roots[0] : roots[1];
                            if (this->verbosity() >= 3)
                                std::cout << Fmt::format("  -- Dogleg step τ={:.2e} (roots: λ_m={:.2e}, λ_p={:.2e})\n", tau, roots[0], roots[1]);

                            // take the Dogleg step by combining the scale Newton and Cauchy point steps
                            Backend::axpy(tau, deltaUU_, deltaU_); // deltaU_ = -(p_U + τ(p_B - p_U))
                            alphaUpdate_ = 1.0;

                            // update the solution (u_k+1 = u_k - αΔu)
                            scaledUpdate_(vars, uLastIter, deltaU_, alphaUpdate_);
                        }
                    }
                    break;
                default:
                    DUNE_THROW(Dune::Exception, "Unknown trust region fallback strategy");
            }
        }

        // full Newton step is within the trust region
        else
        {
            alphaUpdate_ = 1.0;
            deltaU_ = deltaU;

            if (this->verbosity() >= 2)
                std::cout << Fmt::format("-- update inside trust region: ||Δu|| <= Δ ({:.2e} <= {:.2e})\n", newtonStepNorm, delta_);
        }

        // now we have a new proposal step: check if we can accept it or not
        // compute the ratio of the actual reduction to the predicted reduction
        const auto rho = computeTrustRegionRatio_(vars, uLastIter);

        const auto updateNorm = this->linearSolver().norm(deltaU_);

        // update the size of the trust region
        using std::min, std::max;
        if (rho < decrementThreshold_)
        {
            if (this->verbosity() >= 2)
                std::cout << Fmt::format("-- decreasing trust region radius Δ: {:.2e} ==> ", delta_);

            delta_ = max(delta_*decrementFactor_, deltaMin_);

            if (this->verbosity() >= 2)
                std::cout << Fmt::format("{:.2e}\n", delta_);
        }
        else if (rho > incrementThreshold_ && updateNorm == delta_)
        {
            if (this->verbosity() >= 2)
                std::cout << Fmt::format("-- increasing trust region radius Δ: {:.2e} ==> ", delta_);

            delta_ = min(incrementFactor_*delta_, deltaMax_);

            if (this->verbosity() >= 2)
                std::cout << Fmt::format("{:.2e}\n", delta_);
        }
        else
        {
            if (this->verbosity() >= 2)
                std::cout << Fmt::format("-- keeping trust region radius Δ: {:.2e}\n", delta_);
        }

        // accept the proposal step if the actual reduction is above the threshold
        if (rho > eta_)
        {
            if (this->verbosity() >= 2)
                std::cout << Fmt::format("-- good model prediction: ϱ={:.2e} ==> accept: u_k+1 = u_k - Δu\n", rho);

            scaledUpdate_(vars, uLastIter, deltaU_, alphaUpdate_);
            alreadyAssembled_ = false; // reassemble the residual and Jacobian in next iteration
        }

        // reject the proposal step
        else
        {
            if (this->verbosity() >= 2)
                std::cout << Fmt::format("-- ϱ={:.2e} ==> reject: u_k+1 = u_k\n", rho);

            // go to the next iteration with the same solution
            // but updated trust region radius
            alreadyAssembled_ = true;

            this->solutionChanged_(vars, uLastIter);
            this->computeResidualReduction_(vars);
        }
    }

private:
    bool solveLinearSystem_(ResidualVector& deltaU) override
    {
        bool converged = false;
        if (!alreadyAssembled_)
        {
            if (this->verbosity() >= 1)
                std::cout << "-- solve Δu = J^-1 r" << std::endl;

            converged = this->linearSolver().solve(
                this->assembler().jacobian(),
                deltaU,
                this->assembler().residual()
            );

            previousDeltaU_ = deltaU; // store the delta for the next iteration
        }
        else
        {
            if (this->verbosity() >= 1)
                std::cout << "-- set Δu_k+1 = Δu_k" << std::endl;

            deltaU = previousDeltaU_; // set the delta to the previous delta
            converged = true;
        }

        if (this->numSteps_ == 0)
        {
           // initialize the trust region radius
            using std::clamp;
            delta_ = clamp(delta0_*this->initialResidual_, deltaMin_, deltaMax_);
            if (this->verbosity() >= 2)
                std::cout << Fmt::format("-- initialize trust region radius Δ={:.2e} ∈ [{},{}]\n", delta_, deltaMin_, deltaMax_);
        }

        return converged;
    }

    // perform a scaled update of the solution and compute the residual reduction
    void scaledUpdate_(Variables& vars, const SolutionVector& uLastIter, const ResidualVector& deltaU, const Scalar alpha)
    {
        auto uCurrentIter = uLastIter;
        Backend::axpy(-alpha, deltaU, uCurrentIter);
        this->solutionChanged_(vars, uCurrentIter);

        if (enableResidualCriterion_)
            this->computeResidualReduction_(vars);

        // update shift
        this->shift_ = Detail::Newton::maxRelativeShift<Scalar>(uLastIter, uCurrentIter);
        if (this->comm().size() > 1)
            this->shift_ = this->comm().max(this->shift_);
    }

    void maybeAssembleGradientEtc_()
    {
        g_ = 0.0;
        this->assembler().jacobian().mtv(this->assembler().residual(), g_); // g = J^T r
        gNorm_ = this->linearSolver().norm(g_); // ||J^T r||

        w_ = 0.0;
        this->assembler().jacobian().mv(g_, w_); // w = J J^T r
        gTBg_ = w_.dot(w_); // compute g^T B g where B = J^T J and g = J^T r
    }

    void computeCauchyPointUpdate_()
    {
        using std::min;
        const auto tau = gTBg_ < 0.0 ? 1.0 : min(gNorm_*gNorm_*gNorm_/(delta_*gTBg_), 1.0);
        alphaUpdate_ = tau*delta_/gNorm_;
        deltaU_ = g_;

        if (this->verbosity() >= 3)
            std::cout << Fmt::format("  -- Cauchy point update: ||Δu||={:.2e}\n", alphaUpdate_*gNorm_);
    }

    // The pair of roots of a*x^2 + b*x + c = 0 in ascending order
    std::array<Scalar, 2>  quadraticRoots_(const Scalar a, const Scalar b, const Scalar c)
    {
        using std::copysign, std::sqrt;
        const auto temp = -0.5*(b + copysign(1.0, b)*sqrt(b*b - 4.0*a*c));
        const auto sol1 = temp/a;
        const auto sol2 = c/temp;

        using std::min, std::max;
        return {{ min(sol1, sol2), max(sol1, sol2) }};
    }

    Scalar computeTrustRegionRatio_(Variables& vars, const SolutionVector& uLastIter)
    {
        maybeAssembleGradientEtc_();

        // for the actual reduction we need to compute the residual reduction
        auto uCurrentIter = uLastIter;
        Backend::axpy(-alphaUpdate_, deltaU_, uCurrentIter);
        this->solutionChanged_(vars, uCurrentIter);
        this->computeResidualReduction_(vars);

        const auto norm = this->reduction_*this->initialResidual_;
        const auto norm0 = this->lastReduction_*this->initialResidual_;
        const auto actualReduction = 0.5*(norm0*norm0 - norm*norm); // f(u_k) - f(u_k - αΔu) where f = 0.5*||r||^2
        if (this->verbosity() >= 2)
            std::cout << Fmt::format("-- actual reduction f(u)-f(u-Δu)={:.2e}, f(u)={:.2e}, f(u-Δu)={:.2e}\n",
                actualReduction, 0.5*norm0*norm0, 0.5*norm*norm);

        using std::isfinite;
        if (!isfinite(actualReduction))
            return eta_;

        // for the predicted reduction we need to compute reduction predicted by the quadratic model
        // m(0) - m(Δu) = f(u_k) - (f(u_k) + g^T Δu + 0.5*Δu^T B Δu)
        deltaU_ *= alphaUpdate_; // scale the update
        alphaUpdate_ = 1.0; // reset the scaling factor

        this->assembler().jacobian().mv(deltaU_, w_); // w = J Δu
        const auto duTJJTdu = w_.dot(w_); // (Δu J)^T (J Δu) = Δu^T J^T J Δu = Δu^T B Δu
        const auto duTg = deltaU_.dot(g_); // g^T Δu
        const auto predictedReduction = -(-duTg + 0.5*duTJJTdu);
        if (this->verbosity() >= 2)
            std::cout << Fmt::format("-- predicted reduction m(0)-m(Δu)={:.2e}, g^T Δu={:.2e}, Δu^T B Δu={:.2e}\n",
                predictedReduction, duTg, duTJJTdu);

        if (predictedReduction <= 0.0)
            return eta_;

        return actualReduction/predictedReduction;
    }

    bool alreadyAssembled_ = false; //!< flag to check if we need to reassemble the residual and Jacobian
    Scalar delta_; //!< trust region radius
    Scalar delta0_; //!< initial trust region radius
    Scalar deltaMin_; //!< minimum trust region radius
    Scalar deltaMax_; //!< maximum trust region radius
    Scalar eta_; //!< proposal acceptance threshold
    Scalar incrementFactor_; //!< factor to increase the trust region radius
    Scalar decrementFactor_; //!< factor to decrease the trust region radius
    Scalar incrementThreshold_; //!< threshold to increase the trust region radius
    Scalar decrementThreshold_; //!< threshold to decrease the trust region radius
    FallbackStrategy strategy_; //!< trust region fallback strategy

    ResidualVector previousDeltaU_; //!< the update from the previous iteration
    ResidualVector deltaU_; //!< the update for the current iteration
    ResidualVector deltaUU_; //!< temporary vector for the update computation
    Scalar alphaUpdate_; //!< scaling factor for the update step

    Scalar gTBg_; //!< g^T B g where B = J^T J and g = J^T r
    Scalar gNorm_; //!< ||J^T r|| norm of the gradient of the objective function

    ResidualVector g_; //!< gradient of the objective function
    ResidualVector w_; //!< temporary vector for the gradient computation

    bool enableResidualCriterion_ = false; //!< flag to check if we need to compute the residual reduction
};

} // end namespace Dumux

#endif