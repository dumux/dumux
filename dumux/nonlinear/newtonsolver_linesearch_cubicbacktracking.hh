// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Newton solver with line search and cubic backtracking
 */
#ifndef DUMUX_NEWTON_SOLVER_LINESEARCH_CUBIC_BACKTRACKING_HH
#define DUMUX_NEWTON_SOLVER_LINESEARCH_CUBIC_BACKTRACKING_HH

#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \brief An implementation of a Newton solver
 * \tparam Assembler the assembler
 * \tparam LinearSolver the linear solver
 * \tparam Comm the communication object used to communicate with all processes
 * \note If you want to specialize only some methods but are happy with the
 *       defaults of the reference solver, derive your solver from
 *       this class and simply overload the required methods.
 */
template <class Assembler, class LinearSolver,
          class Reassembler = PartialReassembler<Assembler>,
          class Comm = Dune::Communication<Dune::MPIHelper::MPICommunicator> >
class NewtonSolverLinesearchCubicbacktracking
: public NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>
{
    using ParentType = NewtonSolver<Assembler, LinearSolver, Reassembler, Comm>;

protected:
    using Backend = typename ParentType::Backend;
    using SolutionVector = typename ParentType::SolutionVector;
    using ResidualVector = typename ParentType::ResidualVector;

public:
    using typename ParentType::Variables;
    using Communication = Comm;

    NewtonSolverLinesearchCubicbacktracking(
        std::shared_ptr<Assembler> assembler,
        std::shared_ptr<LinearSolver> linearSolver,
        const Communication& comm = Dune::MPIHelper::getCommunication(),
        const std::string& paramGroup = "",
        const std::string& paramGroupName = "Newton",
        int verbosity = 2
    )
    : ParentType(assembler, linearSolver, comm, paramGroup, paramGroupName, verbosity)
    {
        lineSearchMinRelaxationFactor_ = getParam<Scalar>("Newton.LineSearchMinRelaxationFactor", 0.125);
        armijoC1_ = getParam<Scalar>("Newton.LineSearchArmijoC1", 1e-4);
    }

private:
    void lineSearchUpdate_(Variables &vars,
                           const SolutionVector &uLastIter,
                           const ResidualVector &deltaU) override
    {
        const auto& residual = this->assembler().residual();

        Scalar lambda = 1.0;
        auto uCurrentIter = uLastIter;
        Scalar alpha = 1e-4;

        auto ss = residual; ss = 0.0;

        // TODO, in parallel the matrix-vector product requires communication
        if (this->comm().size() > 1)
            DUNE_THROW(Dune::NotImplemented, "Cubic backtracking line search not implemented yet for parallel runs");
        this->assembler().jacobian().mv(deltaU, ss);
        const auto s = residual.dot(ss);

        Backend::axpy(-lambda, deltaU, uCurrentIter);
        this->solutionChanged_(vars, uCurrentIter);
        this->computeResidualReduction_(vars);

        auto norm = this->reduction_*this->initialResidual_;
        const auto norm0 = this->lastReduction_*this->initialResidual_;

        // if the Armijo–Goldstein condition (1st Wolfe condition) is satisfied, we are done
        if (norm*norm <= norm0*norm0 + 2*alpha*lambda*s)
        {
            this->shift_ = Detail::Newton::maxRelativeShift<Scalar>(uLastIter, uCurrentIter);
            if (this->comm().size() > 1)
                this->shift_ = this->comm().max(this->shift_);

            this->endIterMsgStream_ << Fmt::format(", residual reduction {:.4e}->{:.4e}, λ={:.4f}", this->lastReduction_, this->reduction_, lambda);
            return;
        }

        // if not we try to use quadratic backtracking first
        auto mu = lambda;
        auto muNorm = norm;

        uCurrentIter = uLastIter;

        const auto lambda_q = -s*lambda*lambda/(norm*norm - norm0*norm0 - 2*lambda*s);
        lambda = std::clamp(lambda_q, 0.1*lambda, 0.5*lambda);

        if (lambda <= lineSearchMinRelaxationFactor)
            lambda = lineSearchMinRelaxationFactor;

        Backend::axpy(-lambda, deltaU, uCurrentIter);
        this->solutionChanged_(vars, uCurrentIter);
        this->computeResidualReduction_(vars);

        std::cout << Fmt::format("Check quadratic approximation (λ={:.4f})", lambda) << std::endl;

        // if the Armijo–Goldstein condition (1st Wolfe condition) is satisfied, we are done
        norm = this->reduction_*this->initialResidual_;
        if (norm*norm <= norm0*norm0 + 2*alpha*lambda*s || lambda <= lineSearchMinRelaxationFactor)
        {
            this->shift_ = Detail::Newton::maxRelativeShift<Scalar>(uLastIter, uCurrentIter);
            if (this->comm().size() > 1)
                this->shift_ = this->comm().max(this->shift_);

            this->endIterMsgStream_ << Fmt::format(", residual reduction {:.4e}->{:.4e}, λ={:.4f}", this->lastReduction_, this->reduction_, lambda);
            return;
        }

        // if not we try to use cubic backtracking, which is an iterative process
        for (int i = 0; i < 5; ++i)
        {
            uCurrentIter = uLastIter;

            const auto t1 = 0.5*(norm*norm - norm0*norm0) - lambda*s;
            const auto t2 = 0.5*(muNorm*muNorm - norm0*norm0) - mu*s;
            const auto a = (t1/(lambda*lambda) - t2/(mu*mu))/(lambda - mu);
            const auto b = (lambda*t2/(mu*mu) - mu*t1/(lambda*lambda))/(lambda - mu);
            const auto d = std::max(0.0, b*b - 3*a*s);
            const auto lambda_c = a == 0 ? -s/(2*b) : (-b + std::sqrt(d))/(3*a);

            mu = lambda;
            muNorm = norm;

            lambda = std::clamp(lambda_c, 0.1*lambda, 0.5*lambda);

            if (lambda <= lineSearchMinRelaxationFactor)
                lambda = lineSearchMinRelaxationFactor;

            uCurrentIter = uLastIter;

            Backend::axpy(-lambda, deltaU, uCurrentIter);
            this->solutionChanged_(vars, uCurrentIter);
            this->computeResidualReduction_(vars);
            std::cout << Fmt::format("Check cubic approximation (i={}, λ={:.4f}), µ={:.4f})", i, lambda, mu) << std::endl;

            norm = this->reduction_*this->initialResidual_;

            if (norm*norm <= norm0*norm0 + 2*alpha*lambda*s || i == 4 || lambda <= lineSearchMinRelaxationFactor)
            {
                this->shift_ = Detail::Newton::maxRelativeShift<Scalar>(uLastIter, uCurrentIter);
                if (this->comm().size() > 1)
                    this->shift_ = this->comm().max(this->shift_);

                this->endIterMsgStream_ << Fmt::format(", residual reduction {:.4e}->{:.4e}, λ={:.4f}", this->lastReduction_, this->reduction_, lambda);
                return;
            }
        }
    }

    Scalar lineSearchMinRelaxationFactor_;
    Scalar armijoC1_;
};

} // end namespace Dumux

#endif
