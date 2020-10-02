// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Linear
 * \brief Dumux preconditioners for iterative solvers
 */
#ifndef DUMUX_LINEAR_PRECONDITIONERS_HH
#define DUMUX_LINEAR_PRECONDITIONERS_HH

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/indices.hh>
#include <dune/common/version.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/paamg/amg.hh>

#if HAVE_UMFPACK
#include <dune/istl/umfpack.hh>
#endif

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/matrix.hh>
#include <dumux/linear/istlsolverregistry.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A preconditioner based on the Uzawa algorithm for saddle-point problems of the form
 * \f$
 \begin{pmatrix}
    A & B \\
    C & D
 \end{pmatrix}

 \begin{pmatrix}
    u\\
    p
 \end{pmatrix}

 =

 \begin{pmatrix}
    f\\
    g
 \end{pmatrix}
  * \f$
 *
 * This preconditioner is especially suited for solving the incompressible (Navier-)Stokes equations.
 * Here, \f$D = 0\f$ and \f$B = C^T\f$ if \f$\rho = 1\f$.
 * We do not expect good convergence if energy or mass transport is considered.
 *
 * See: Benzi, M., Golub, G. H., & Liesen, J. (2005). Numerical solution of saddle point problems. Acta numerica, 14, 1-137 \cite benzi2005 and <BR>
 *      Ho, N., Olson, S. D., & Walker, H. F. (2017). Accelerating the Uzawa algorithm. SIAM Journal on Scientific Computing, 39(5), S461-S476 \cite ho2017
 *
 * \tparam M Type of the matrix.
 * \tparam X Type of the update.
 * \tparam Y Type of the defect.
 * \tparam l Preconditioner block level (for compatibility reasons, unused).
 */
template<class M, class X, class Y, int l = 1>
class SeqUzawa : public Dune::Preconditioner<X,Y>
{
    static_assert(Dumux::isMultiTypeBlockMatrix<M>::value && M::M() == 2 && M::N() == 2, "SeqUzawa expects a 2x2 MultiTypeBlockMatrix.");
    static_assert(l== 1, "SeqUzawa expects a block level of 1.");

    using A = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_0][Dune::Indices::_0])>;
    using U = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_0])>;

    using Comm = Dune::Amg::SequentialInformation;
    using LinearOperator = Dune::MatrixAdapter<A, U, U>;
    using Smoother = Dune::SeqSSOR<A, U, U>;
    using AMGSolverForA = Dune::Amg::AMG<LinearOperator, U, Smoother, Comm>;

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = M;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;
    //! \brief Scalar type underlying the field_type.
    using scalar_field_type = Dune::Simd::Scalar<field_type>;

    /*!
     * \brief Constructor
     *
     * \param mat The matrix to operate on.
     * \param params Collection of paramters.
     */
#if DUNE_VERSION_GT(DUNE_ISTL,2,7)
    SeqUzawa(const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& op, const Dune::ParameterTree& params)
    : matrix_(op->getmat())
#else
    SeqUzawa(const M& mat, const Dune::ParameterTree& params)
    : matrix_(mat)
#endif
    , numIterations_(params.get<std::size_t>("iterations"))
    , relaxationFactor_(params.get<scalar_field_type>("relaxation"))
    , verbosity_(params.get<int>("verbosity"))
    , paramGroup_(params.get<std::string>("ParameterGroup"))
    , useDirectVelocitySolverForA_(getParamFromGroup<bool>(paramGroup_, "LinearSolver.Preconditioner.DirectSolverForA", false))
    {
        const bool determineRelaxationFactor = getParamFromGroup<bool>(paramGroup_, "LinearSolver.Preconditioner.DetermineRelaxationFactor", true);

        // AMG is needed for determination of omega
        if (determineRelaxationFactor || !useDirectVelocitySolverForA_)
            initAMG_(params);

        if (useDirectVelocitySolverForA_)
            initUMFPack_();

        if (determineRelaxationFactor)
            relaxationFactor_ = estimateOmega_();
    }

    /*!
     * \brief Prepare the preconditioner.
     */
    virtual void pre(X& x, Y& b) {}

    /*!
     * \brief Apply the preconditioner
     *
     * \param update The update to be computed.
     * \param currentDefect The current defect.
     */
    virtual void apply(X& update, const Y& currentDefect)
    {
        using namespace Dune::Indices;

        auto& A = matrix_[_0][_0];
        auto& B = matrix_[_0][_1];
        auto& C = matrix_[_1][_0];
        auto& D = matrix_[_1][_1];

        const auto& f = currentDefect[_0];
        const auto& g = currentDefect[_1];
        auto& u = update[_0];
        auto& p = update[_1];

        // incorporate Dirichlet cell values
        // TODO: pass Dirichlet constraint handler from outside
        for (std::size_t i = 0; i < D.N(); ++i)
        {
            const auto& block = D[i][i];
            for (auto rowIt = block.begin(); rowIt != block.end(); ++rowIt)
                if (Dune::FloatCmp::eq<scalar_field_type>(rowIt->one_norm(), 1.0))
                    p[i][rowIt.index()] = g[i][rowIt.index()];
        }

        // the actual Uzawa iteration
        for (std::size_t k = 0; k < numIterations_; ++k)
        {
            // u_k+1 = u_k + Q_A^−1*(f − (A*u_k + B*p_k)),
            auto uRhs = f;
            A.mmv(u, uRhs);
            B.mmv(p, uRhs);
            auto uIncrement = u;
            applySolverForA_(uIncrement, uRhs);
            u += uIncrement;

            // p_k+1 = p_k + omega*(g - C*u_k+1 - D*p_k)
            auto pIncrement = g;
            C.mmv(u, pIncrement);
            D.mmv(p, pIncrement);
            pIncrement *= relaxationFactor_;
            p += pIncrement;

            if (verbosity_ > 1)
            {
                std::cout << "Uzawa iteration " << k
                << ", residual: " << uRhs.two_norm() + pIncrement.two_norm()/relaxationFactor_ << std::endl;
            }
        }
    }

    /*!
     * \brief Clean up.
     */
    virtual void post(X& x) {}

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
    }

private:

    void initAMG_(const Dune::ParameterTree& params)
    {
        using namespace Dune::Indices;
        auto linearOperator = std::make_shared<LinearOperator>(matrix_[_0][_0]);
        amgSolverForA_ = std::make_unique<AMGSolverForA>(linearOperator, params);
    }

    void initUMFPack_()
    {
#if HAVE_UMFPACK
            using namespace Dune::Indices;
            umfPackSolverForA_ = std::make_unique<Dune::UMFPack<A>>(matrix_[_0][_0]);
#else
            DUNE_THROW(Dune::InvalidStateException, "UMFPack not available. Use LinearSolver.Preconditioner.DirectVelocitySolver = false.");
#endif
    }

    /*!
     * \brief Estimate the relaxation factor omega
     *
     * The optimal relaxation factor is omega = 2/(lambdaMin + lambdaMax), where lambdaMin and lambdaMax are the smallest and largest
     * eigenvalues of the Schur complement -C*Ainv*B (assuming D = 0).
     * lambdaMax can be easily determined using the power iteration algorithm (https://en.wikipedia.org/wiki/Power_iteration) and lambdaMin
     * could be estimated in a similar manner. We do not consider lambdaMin because for certain cases, e.g., when C contains some rows of zeroes only,
     * this estimate will fail.
     *
     * Instead we assume that lambdaMin is sufficiently close to lambdaMax such that omega = 1/lambdaMax.
     * This seems to work rather well for various applications.
     * We will underestimate omega by a factor of 2 in the worst case (i.e, lambdaMin = 0).
     *
     * When facing convergence issues, you may set LinearSolver.Preconditioner.Verbosity = 1 to see the estimate of lambdaMax.
     * In a new simulation run, you can then set LinearSolver.Preconditioner.DetermineRelaxationFactor = false and set some other value
     * for LinearSolver.Preconditioner.Relaxation based on the estimate of lambdaMax.
     *
     * See: Benzi, M., Golub, G. H., & Liesen, J. (2005). Numerical solution of saddle point problems. Acta numerica, 14, 1-137.
     */
    scalar_field_type estimateOmega_() const
    {
        using namespace Dune::Indices;
        auto& A = matrix_[_0][_0];
        auto& B = matrix_[_0][_1];
        auto& C = matrix_[_1][_0];

        U x(A.M());
        x = 1.0;

        scalar_field_type omega = 0.0;
        scalar_field_type lambdaMax = 0.0;

        static const auto iterations = Dumux::getParamFromGroup<std::size_t>(paramGroup_, "LinearSolver.Preconditioner.PowerLawIterations", 5);

        // apply power iteration x_k+1 = M*x_k/|M*x_k| for the matrix M = -C*Ainv*B
        for (std::size_t i = 0; i < iterations; ++i)
        {
            // bx = B*x
            U bx(x.size());
            B.mv(x, bx);

            // ainvbx = Ainv*(B*x)
            auto ainvbx = x;
            applySolverForA_(ainvbx, bx);

            // v = M*x = -C*(Ainv*B*x)
            U v(x.size());
            C.mv(ainvbx, v);
            v *= -1.0;

            // eigenvalue lambdaMax = xt*M*x/(xt*x) = xt*v/(xt*x);
            lambdaMax = x.dot(v)/(x.dot(x));

            // relaxation factor omega = 1/lambda;
            omega = 1.0/lambdaMax;

            // new iterate x = M*x/|M*x| = v/|v|
            x = v;
            x /= v.two_norm();
        }

        if (verbosity_ > 0)
        {
            std::cout << "\n*** Uzawa Preconditioner ***" << std::endl;
            std::cout << "Estimating relaxation factor based on Schur complement" << std::endl;
            std::cout << "Largest estimated eigenvalue lambdaMax = " << lambdaMax << std::endl;
            std::cout << "Relaxation factor omega = 1/lambdaMax = " << omega << std::endl;
        }

        return omega;
    }

    template<class Sol, class Rhs>
    void applySolverForA_(Sol& sol, Rhs& rhs) const
    {
        if (useDirectVelocitySolverForA_)
        {
#if HAVE_UMFPACK
            Dune::InverseOperatorResult res;
            umfPackSolverForA_->apply(sol, rhs, res);
#endif
        }
        else
        {
            amgSolverForA_->pre(sol, rhs);
            amgSolverForA_->apply(sol, rhs);
            amgSolverForA_->post(sol);
        }
    }

    //! \brief The matrix we operate on.
    const M& matrix_;
    //! \brief The number of steps to do in apply
    const std::size_t numIterations_;
    //! \brief The relaxation factor to use
    scalar_field_type relaxationFactor_;
    //! \brief The verbosity level
    const int verbosity_;

    std::unique_ptr<AMGSolverForA> amgSolverForA_;
#if HAVE_UMFPACK
    std::unique_ptr<Dune::UMFPack<A>> umfPackSolverForA_;
#endif
    const std::string paramGroup_;
    const bool useDirectVelocitySolverForA_;
};

DUMUX_REGISTER_PRECONDITIONER("uzawa", Dumux::MultiTypeBlockMatrixPreconditionerTag, Dune::defaultPreconditionerBlockLevelCreator<Dumux::SeqUzawa, 1>());

} // end namespace Dumux

#endif
