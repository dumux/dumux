// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Newton
 * \brief Newton solver with a line search by quadratic and cubic backtracking
 */
#ifndef DUMUX_NEWTON_SOLVER_LINESEARCH_CUBIC_BACKTRACKING_HH
#define DUMUX_NEWTON_SOLVER_LINESEARCH_CUBIC_BACKTRACKING_HH

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

#include <dumux/common/parameters.hh>
#include <dumux/io/format.hh>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux {

/*!
 * \ingroup Newton
 * \brief A Newton solver with a line search by quadratic and cubic backtracking
 *
 * The update \f$ \textbf{u}_{k+1} = \textbf{u}_k - \alpha \Delta \textbf{u} \f$ is accepted if
 * \f$ g(\alpha) = \frac{1}{2} \Vert F(\textbf{u}_k - \alpha \Delta \textbf{u}) \Vert_2^2 \f$
 * satisfies the sufficient decrease (Armijo) condition \f$ g(\alpha) \leq g(0) + c_1 \alpha g^\prime(0) \f$.
 * Starting from \f$ \alpha = 1 \f$, a rejected step size is replaced by the minimizer of the quadratic
 * (first rejection) or cubic (subsequent rejections) polynomial interpolating \f$ g \f$,
 * restricted to \f$ [0.1\alpha, 0.5\alpha] \f$, see \cite NocedalWright2006 (Section 3.5).
 * If the residual is not finite for a trial step size, the step size is halved.
 * A step size that reaches the lower bound `Newton.LineSearchMinRelaxationFactor` (default \f$ 0.125 \f$)
 * is accepted even if it does not satisfy the sufficient decrease condition.
 * The constant \f$ c_1 \f$ is set by `Newton.LineSearchArmijoC1` (default \f$ 10^{-4} \f$).
 * The line search is always enabled, independent of `Newton.UseLineSearch`.
 *
 * \tparam Assembler the assembler
 * \tparam LinearSolver the linear solver
 * \tparam Reassembler the reassembler
 * \tparam Comm the communication object used to communicate with all processes
 */
template <class Assembler, class LinearSolver,
          class Reassembler = PartialReassembler<Assembler>,
          class Comm = Dune::Communication<Dune::MPIHelper::MPICommunicator> >
class NewtonSolverLinesearchCubicbacktracking
: public NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>
{
    using ParentType = NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>;
    using Scalar = typename Assembler::Scalar;

protected:
    using Backend = typename ParentType::Backend;
    using SolutionVector = typename ParentType::SolutionVector;
    using ResidualVector = typename ParentType::ResidualVector;

public:
    using typename ParentType::Variables;
    using Communication = Comm;

    NewtonSolverLinesearchCubicbacktracking(std::shared_ptr<Assembler> assembler,
                                            std::shared_ptr<LinearSolver> linearSolver,
                                            const Communication& comm = Dune::MPIHelper::getCommunication(),
                                            const std::string& paramGroup = "",
                                            const std::string& paramGroupName = "Newton",
                                            int verbosity = 2)
    : ParentType(assembler, linearSolver, comm, paramGroup, paramGroupName, verbosity)
    {
        this->setUseLineSearch(true);
        minRelaxationFactor_ = getParamFromGroup<Scalar>(paramGroup, paramGroupName + ".LineSearchMinRelaxationFactor", 0.125);
        armijoC1_ = getParamFromGroup<Scalar>(paramGroup, paramGroupName + ".LineSearchArmijoC1", 1e-4);
    }

private:
    void lineSearchUpdate_(Variables& vars,
                           const SolutionVector& uLastIter,
                           const ResidualVector& deltaU) override
    {
        using std::isfinite;

        const Scalar residualNorm0 = this->lastReduction_*this->initialResidual_;
        const Scalar g0 = 0.5*residualNorm0*residualNorm0;
        // g'(0) = -r^T J deltaU = -|r|^2 because the Newton update solves J deltaU = r
        const Scalar slope = -residualNorm0*residualNorm0;

        Scalar lambda = 1.0;
        Scalar lambdaPrev = 0.0;
        Scalar gPrev = 0.0;
        bool hasPrevTrial = false;
        auto uCurrentIter = uLastIter;
        while (true)
        {
            Backend::axpy(-lambda, deltaU, uCurrentIter);
            this->solutionChanged_(vars, uCurrentIter);
            this->computeResidualReduction_(vars);

            const Scalar g = 0.5*this->residualNorm_*this->residualNorm_;
            if (g <= g0 + armijoC1_*lambda*slope || lambda <= minRelaxationFactor_)
            {
                this->endIterMsgStream_ << Fmt::format(", residual reduction {:.4e}->{:.4e}@lambda={:.4f}", this->lastReduction_, this->reduction_, lambda);
                return;
            }

            Scalar lambdaNext = 0.5*lambda;
            if (isfinite(g))
            {
                const Scalar lambdaModel = hasPrevTrial ? cubicMinimizer_(g0, slope, lambda, g, lambdaPrev, gPrev)
                                                        : quadraticMinimizer_(g0, slope, lambda, g);
                if (isfinite(lambdaModel))
                    lambdaNext = lambdaModel;

                lambdaPrev = lambda;
                gPrev = g;
                hasPrevTrial = true;
            }

            lambda = std::max<Scalar>(std::clamp<Scalar>(lambdaNext, 0.1*lambda, 0.5*lambda), minRelaxationFactor_);
            uCurrentIter = uLastIter;
        }
    }

    // minimizer of the quadratic q with q(0) = g0, q'(0) = slope, q(lambda) = g
    static Scalar quadraticMinimizer_(Scalar g0, Scalar slope, Scalar lambda, Scalar g)
    { return -slope*lambda*lambda/(2.0*(g - g0 - slope*lambda)); }

    // minimizer of the cubic c with c(0) = g0, c'(0) = slope, c(lambda) = g, c(lambdaPrev) = gPrev
    static Scalar cubicMinimizer_(Scalar g0, Scalar slope, Scalar lambda, Scalar g, Scalar lambdaPrev, Scalar gPrev)
    {
        using std::sqrt;
        const Scalar t = (g - g0 - slope*lambda)/(lambda*lambda);
        const Scalar tPrev = (gPrev - g0 - slope*lambdaPrev)/(lambdaPrev*lambdaPrev);
        const Scalar a = (t - tPrev)/(lambda - lambdaPrev);
        const Scalar b = (lambda*tPrev - lambdaPrev*t)/(lambda - lambdaPrev);
        if (a == 0.0)
            return -slope/(2.0*b);

        const Scalar discriminant = b*b - 3.0*a*slope;
        if (discriminant < 0.0)
            return 0.5*lambda;
        else if (b <= 0.0)
            return (-b + sqrt(discriminant))/(3.0*a);
        else // same root as above without cancellation
            return -slope/(b + sqrt(discriminant));
    }

    Scalar minRelaxationFactor_;
    Scalar armijoC1_;
};

} // end namespace Dumux

#endif
