// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Nonlinear
 * \brief Levenberg-Marquardt algorithm for solving nonlinear least squares problems
 */
#ifndef DUMUX_NONLINEAR_LEASTSQUARES_HH
#define DUMUX_NONLINEAR_LEASTSQUARES_HH

#include <iostream>
#include <cmath>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrix.hh>
#if HAVE_SUITESPARSE_CHOLMOD
#include <dune/istl/cholmod.hh>
#endif // HAVE_SUITESPARSE_CHOLMOD
#include <dune/istl/io.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/io/format.hh>
#include <dumux/nonlinear/newtonsolver.hh>

namespace Dumux::Optimization {

//! a solver base class
template<class Variables>
class Solver
{
public:
    virtual ~Solver() = default;
    virtual bool apply(Variables& vars) = 0;
};

} // end namespace Dumux::Optimization

#ifndef DOXYGEN // exclude from doxygen
// helper code for curve fitting
namespace Dumux::Optimization::Detail {

template<class T, class F>
class Assembler
{
public:
    using Scalar = T;
    using JacobianMatrix = Dune::Matrix<T>;
    using ResidualType = Dune::BlockVector<T>;
    using SolutionVector = Dune::BlockVector<T>;
    using Variables = Dune::BlockVector<T>;

    /*!
     * \brief Assembler for the Levenberg-Marquardt linear system: [J^T J + λ diag(J^T J)] x = J^T r
     * \param f The residual function r = f(x) to minimize the norm of (\f$ f: \mathbb{R}^n \mapsto \mathbb{R}^m \f$)
     * \param x0 Initial guess for the parameters (vector of size \f$n\f$)
     * \param residualSize The number of residuals (\f$m\f$)
     */
    Assembler(const F& f, const SolutionVector& x0, std::size_t residualSize)
    : prevSol_(x0), f_(f), solSize_(x0.size()), residualSize_(residualSize)
    , JT_(solSize_, residualSize_), regularizedJTJ_(solSize_, solSize_)
    , residual_(residualSize_), projectedResidual_(solSize_)
    {
        printMatrix_ = getParamFromGroup<bool>("LevenbergMarquardt", "PrintMatrix", false);
        baseEpsilon_ = getParamFromGroup<Scalar>("LevenbergMarquardt", "BaseEpsilon", 1e-1);
    }

    //! initialize the linear system variables
    void setLinearSystem()
    {
        std::cout << "Setting up linear system with " << solSize_ << " variables and "
                  << residualSize_ << " residuals." << std::endl;

        JT_ = 0.0;
        regularizedJTJ_ = 0.0;
        residual_ = 0.0;
        projectedResidual_ = 0.0;
    }

    //! interface for the Newton scheme: this is J^T J + λ diag(J^T J)
    JacobianMatrix& jacobian() { return regularizedJTJ_; }

    //! interface for the Newton scheme: this is J^T r
    ResidualType& residual() { return projectedResidual_; }

    //! regularization parameter λ
    void setLambda(const T lambda) { lambda_ = lambda; }
    T lambda() { return lambda_; }

    //! assemble the right hand side of the linear system
    void assembleResidual(const SolutionVector& x)
    {
        projectedResidual_ = 0.0;
        JT_ = 0.0;

        // assemble residual, J^T, and projected residual, J^T r
        for (auto rowIt = JT_.begin(); rowIt != JT_.end(); ++rowIt)
        {
            const auto paramIdx = rowIt.index();
            const auto residual = f_(x);
            auto p = x;
            const auto eps = baseEpsilon_;
            p[paramIdx] = x[paramIdx] + eps;
            auto deflectedResidual = f_(p);
            deflectedResidual -= residual;
            deflectedResidual /= eps;
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                *colIt += deflectedResidual[colIt.index()];

            projectedResidual_[paramIdx] = residual * deflectedResidual;
        }

        if (printMatrix_)
        {
            std::cout << std::endl << "J^T = " << std::endl;
            Dune::printmatrix(std::cout, JT_, "", "");
        }
    }

    //! assemble the left hand side of the linear system
    void assembleJacobianAndResidual(const SolutionVector& x)
    {
        assembleResidual(x);

        regularizedJTJ_ = 0.0;
        for (auto rowIt = regularizedJTJ_.begin(); rowIt != regularizedJTJ_.end(); ++rowIt)
        {
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
            {
                for (int j = 0; j < residualSize_; ++j)
                    *colIt += JT_[rowIt.index()][j]*JT_[colIt.index()][j];

                if (rowIt.index() == colIt.index())
                {
                    *colIt += lambda_*(*colIt);

                    if (*colIt == 0.0)
                        *colIt = 1e-6*lambda_;
                }
            }
        }

        if (printMatrix_)
        {
            std::cout << std::endl << "J^T J + λ diag(J^T J) = " << std::endl;
            Dune::printmatrix(std::cout, regularizedJTJ_, "", "");
        }
    }

    //! this is the actual cost unction we are minimizing MSE = ||f_(x)||^2
    T computeMeanSquaredError(const SolutionVector& x) const
    { return f_(x).two_norm2(); }

    //! initial guess to be able to reset the solver
    const SolutionVector& prevSol() const
    { return prevSol_; }

private:
    SolutionVector prevSol_; // initial guess
    const F& f_; //!< the residual function to minimize the norm of
    std::size_t solSize_, residualSize_;

    JacobianMatrix JT_; // J^T
    JacobianMatrix regularizedJTJ_; // J^T J + λ diag(J^T J)
    ResidualType residual_; // r = f_(x)
    ResidualType projectedResidual_; // J^T r

    Scalar lambda_ = 0.0; // regularization parameter
    Scalar baseEpsilon_; // for numerical differentiation

    bool printMatrix_ = false;
};

#if HAVE_SUITESPARSE_CHOLMOD
template<class Matrix, class Vector>
class CholmodLinearSolver
{
public:
    void setResidualReduction(double residualReduction) {}

    bool solve(const Matrix& A, Vector& x, const Vector& b) const
    {
        Dune::Cholmod<Vector> solver; // matrix is symmetric
        solver.setMatrix(A);
        Dune::InverseOperatorResult r;
        auto bCopy = b;
        auto xCopy = x;
        solver.apply(xCopy, bCopy, r);
        checkResult_(xCopy, r);
        if (!r.converged )
            DUNE_THROW(NumericalProblem, "Linear solver did not converge.");
        x = xCopy;
        return r.converged;
    }

    auto norm(const Vector& residual) const
    { return residual.two_norm(); }

private:
    void checkResult_(Vector& x, Dune::InverseOperatorResult& result) const
    {
        flatVectorForEach(x, [&](auto&& entry, std::size_t){
            using std::isnan, std::isinf;
            if (isnan(entry) || isinf(entry))
                result.converged = false;
        });
    }
};
#endif // HAVE_SUITESPARSE_CHOLMOD

} // end namespace Dumux::Optimization::Detail
#endif // DOXYGEN

namespace Dumux::Optimization::Detail {

/*!
 * \ingroup Nonlinear
 * \brief A nonlinear least squares solver with \f$n\f$ model parameters and \f$m\f$ observations
 * \tparam Assembler a class that assembles the linear system of equations in each iteration
 * \tparam LinearSolver a class that solves the linear system of equations
 * \note The assembler has to assemble the actual least squares problem (i.e. the normal equations)
 */
template<class Assembler, class LinearSolver>
class LevenbergMarquardt : private Dumux::NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler, Dune::Communication<Dune::No_Comm>>
{
    using ParentType = Dumux::NewtonSolver<Assembler, LinearSolver, DefaultPartialReassembler, Dune::Communication<Dune::No_Comm>>;
    using Scalar = typename Assembler::Scalar;
    using Variables = typename Assembler::Variables;
    using SolutionVector = typename Assembler::SolutionVector;
    using ResidualVector = typename Assembler::ResidualType;
    using Backend = typename ParentType::Backend;
public:
    using ParentType::ParentType;

    LevenbergMarquardt(std::shared_ptr<Assembler> assembler,
                       std::shared_ptr<LinearSolver> linearSolver,
                       int verbosity = 0)
    : ParentType(assembler, linearSolver, {}, "", "LevenbergMarquardt", verbosity)
    {
        assembler->setLambda(getParamFromGroup<Scalar>(ParentType::paramGroup(), "LevenbergMarquardt.LambdaInitial", 0.001));
        maxLambdaDivisions_ = getParamFromGroup<std::size_t>(ParentType::paramGroup(), "LevenbergMarquardt.MaxLambdaDivisions", 10);

        this->setUseLineSearch();
        this->setMaxRelativeShift(1e-3);
    }

    /*!
     * \brief Solve the nonlinear least squares problem
     * \param vars instance of the `Variables` class representing a numerical
     *             solution, defining primary and possibly secondary variables
     *             and information on the time level.
     * \return bool true if the solver converged
     * \post If converged, the given `Variables` will represent the solution. If the solver
     *       does not converge, the `Variables` will represent the best solution so far (lowest value of the objective function)
     */
    bool apply(Variables& vars) override
    {
        // try solving the non-linear system
        for (std::size_t i = 0; i <= maxLambdaDivisions_; ++i)
        {
            // linearize & solve
            const bool converged = ParentType::apply(vars);

            if (converged)
                return true;

            else if (!converged && i < maxLambdaDivisions_)
            {
                if (this->verbosity() >= 1)
                    std::cout << Fmt::format("LevenbergMarquardt solver did not converge with λ = {:.2e}. ", ParentType::assembler().lambda())
                                << Fmt::format("Retrying with λ = {:.2e}.\n", ParentType::assembler().lambda() * 9);

                ParentType::assembler().setLambda(ParentType::assembler().lambda() * 9);
            }

            // convergence criterium not fulfilled
            // return best solution so far
            else
            {
                this->solutionChanged_(vars, bestSol_);
                if (this->verbosity() >= 1)
                    std::cout << Fmt::format("Choose best solution so far with a MSE of {:.4e}", minResidual_) << std::endl;

                return false;
            }
        }

        return false;
    }

private:
    void newtonEndStep(Variables &vars, const SolutionVector &uLastIter) override
    {
        ParentType::newtonEndStep(vars, uLastIter);

        if (ParentType::reduction_ > ParentType::lastReduction_)
        {
            if (ParentType::assembler().lambda() < 1e5)
            {
                ParentType::assembler().setLambda(ParentType::assembler().lambda() * 9);
                if (this->verbosity() > 0)
                    std::cout << "-- Increasing λ to " << ParentType::assembler().lambda();
            }
            else
            {
                if (this->verbosity() > 0)
                    std::cout << "-- Keeping λ at " << ParentType::assembler().lambda();
            }
        }
        else
        {
            if (ParentType::reduction_ < 0.1 && ParentType::assembler().lambda() > 1e-5)
            {
                ParentType::assembler().setLambda(ParentType::assembler().lambda() / 11.0);
                if (this->verbosity() > 0)
                    std::cout << "-- Decreasing λ to " << ParentType::assembler().lambda();
            }
            else
            {
                if (this->verbosity() > 0)
                    std::cout << "-- Keeping λ at " << ParentType::assembler().lambda();
            }
        }

        if (this->verbosity() > 0)
            std::cout << ", error reduction: " << ParentType::reduction_
                      << " and mean squared error: " <<  ParentType::residualNorm_ << std::endl;

        // store best solution
        if (ParentType::residualNorm_ < minResidual_)
        {
            minResidual_ = ParentType::residualNorm_;
            bestSol_ = uLastIter;
        }
    }

    void lineSearchUpdate_(Variables& vars,
                           const SolutionVector& uLastIter,
                           const ResidualVector& deltaU) override
    {
        alpha_ = 1.0;
        auto uCurrentIter = uLastIter;

        while (true)
        {
            Backend::axpy(-alpha_, deltaU, uCurrentIter);
            this->solutionChanged_(vars, uCurrentIter);

            ParentType::residualNorm_ = ParentType::assembler().computeMeanSquaredError(Backend::dofs(vars));
            if (ParentType::numSteps_ == 0)
                ParentType::initialResidual_ = ParentType::assembler().computeMeanSquaredError(ParentType::assembler().prevSol());
            ParentType::reduction_ = ParentType::residualNorm_ / ParentType::initialResidual_;

            if (ParentType::reduction_ < ParentType::lastReduction_ || alpha_ <= 0.001)
            {
                ParentType::endIterMsgStream_ << Fmt::format(", residual reduction {:.4e}->{:.4e}@α={:.4f}", ParentType::lastReduction_, ParentType::reduction_, alpha_);
                return;
            }

            // try with a smaller update and reset solution
            alpha_ *= 0.5;
            uCurrentIter = uLastIter;
        }
    }

    Scalar minResidual_ = std::numeric_limits<Scalar>::max();
    SolutionVector bestSol_;
    std::size_t maxLambdaDivisions_ = 10;
    Scalar alpha_ = 1.0;
};

#if HAVE_SUITESPARSE_CHOLMOD
template<class T, class F>
class NonlinearLeastSquaresSolver : public Solver<typename Optimization::Detail::Assembler<T, F>::Variables>
{
    using Assembler = Optimization::Detail::Assembler<T, F>;
    using LinearSolver = Optimization::Detail::CholmodLinearSolver<typename Assembler::JacobianMatrix, typename Assembler::ResidualType>;
    using Optimizer = Optimization::Detail::LevenbergMarquardt<Assembler, LinearSolver>;
public:
    using Variables = typename Assembler::Variables;

    NonlinearLeastSquaresSolver(const F& f, const Dune::BlockVector<T>& x0, std::size_t size)
    : solver_(std::make_unique<Optimizer>(std::make_shared<Assembler>(f, x0, size), std::make_shared<LinearSolver>(), 2))
    {}

    bool apply(Variables& vars) override
    { return solver_->apply(vars); }

private:
    std::unique_ptr<Optimizer> solver_;
};
#endif // HAVE_SUITESPARSE_CHOLMOD

} // end namespace Dumux::Optimization::Detail

namespace Dumux::Optimization {

#if HAVE_SUITESPARSE_CHOLMOD

/*!
 * \ingroup Nonlinear
 * \brief Creates a nonlinear least squares solver with \f$n\f$ model parameters and \f$m\f$ observations
 * \param f The residual functional to minimize the norm of (\f$ f: \mathbb{R}^n \mapsto \mathbb{R}^m \f$)
 * \param x0 Initial guess for the parameters (vector of size \f$n\f$)
 * \param size The number of observations (\f$m\f$)
 * \return a unique pointer to the nonlinear least squares solver with a method `apply(Variables& vars)`
 * \note The solver can be configured through the parameter group `LevenbergMarquardt`
 */
template<class T, class F>
auto makeNonlinearLeastSquaresSolver(const F& f, const Dune::BlockVector<T>& x0, std::size_t size)
-> std::unique_ptr<Solver<Dune::BlockVector<T>>>
{ return std::make_unique<Optimization::Detail::NonlinearLeastSquaresSolver<T, F>>(f, x0, size); }

#endif // HAVE_SUITESPARSE_CHOLMOD

} // end namespace Dumux::Optimization

#endif
