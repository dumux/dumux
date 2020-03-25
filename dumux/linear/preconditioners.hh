// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \brief Dumux precondioners for iterative solvers
 */
#ifndef DUMUX_PRECONDITIONERS_HH
#define DUMUX_PRECONDITIONERS_HH

#include <dune/common/indices.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/umfpack.hh>
#include <dune/istl/common/registry.hh>
#include <dune/istl/solverregistry.hh>

#include <dune/common/exceptions.hh>
#include <dumux/common/parameters.hh>

#define DUMUX_REGISTER_PRECONDITIONER(name, ...)                \
  DUNE_REGISTRY_PUT(DumuxPreconditionerTag, name, __VA_ARGS__)

struct DumuxPreconditionerTag {};


namespace Dune {

template<class M, class X, class Y, int l = 1>
class NewSeqUzawa : public Dune::Preconditioner<X,Y>
{
    using VelocityMatrix = typename std::remove_reference<decltype(std::declval<M>()[Dune::Indices::_0][Dune::Indices::_0])>::type;
    using VelocityVector = typename std::remove_reference<decltype(std::declval<X>()[Dune::Indices::_0])>::type;

    using Comm = Dune::Amg::SequentialInformation;
    using LinearOperator = Dune::MatrixAdapter<VelocityMatrix, VelocityVector, VelocityVector>;
    using ScalarProduct = Dune::SeqScalarProduct<VelocityVector>;
    using Smoother = Dune::SeqSSOR<VelocityMatrix, VelocityVector, VelocityVector>;
    using VelocityAMG = Dune::Amg::AMG<LinearOperator, VelocityVector, Smoother, Comm>;
    using VelocityUMFPack = Dune::UMFPack<VelocityMatrix>;

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = M;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;
    //! \brief scalar type underlying the field_type
    using scalar_field_type = Dune::Simd::Scalar<field_type>;

    /*! \brief Constructor.
     *
     *   constructor gets all parameters to operate the prec.
     *   \param mat The matrix to operate on.
     *   \param n The number of iterations to perform.
     *   \param w The relaxation factor.
     */
    NewSeqUzawa (const std::shared_ptr<const AssembledLinearOperator<M,X,Y>>& mat, const ParameterTree& params)
    : A_(mat->getmat())
    , n_(params.get<std::size_t>("iterations"))
    , w_(params.get<scalar_field_type>("relaxation"))
    , paramGroup_(params.get<std::string>("ParameterGroup"))
    , inexact_(Dumux::getParamFromGroup<bool>(paramGroup_, "LinearSolver.InexactVelocitySolver", false))
    {
        init_(params);
    }

    /*!
     *   \brief Prepare the preconditioner.
     *
     *   \copydoc Preconditioner::pre(X&,Y&)
     */
    virtual void pre (X& x, Y& b) {}

    /*!
     *   \brief Apply the preconditioner
     *
     *   \copydoc Preconditioner::apply(X&,const Y&)
     */
    virtual void apply (X& v, const Y& d)
    {
        using namespace Dune::Indices;

        auto& A = A_[_0][_0];
        auto& Bt = A_[_0][_1];
        auto& B = A_[_1][_0];
        auto& C = A_[_1][_1];

        const auto& f = d[_0];
        const auto& g = d[_1];
        auto& velocity = v[_0];
        auto& pressure = v[_1];
        for (int i = 0; i < C.N(); ++i)
        {
            if (std::abs(C[i][i].frobenius_norm() - 1.0) < 1e-14)
                pressure[i] = g[i];
        }

        for (int k = 0; k < n_; ++k)
        {
            // u_k+1 = u_k + Q_A^−1*(f − (A*u_k + Bt*p_k)),
            auto vrhs = f;
            A.mmv(velocity, vrhs);
            Bt.mmv(pressure, vrhs);
            auto vUpdate = velocity;
            if (inexact_)
            {
                velocityAMG_->pre(vUpdate, vrhs);
                velocityAMG_->apply(vUpdate, vrhs);
                velocityAMG_->post(vUpdate);
            }
            else
            {
                Dune::InverseOperatorResult res;
                velocityUMFPack_->apply(vUpdate, vrhs, res);
            }
            velocity += vUpdate;

            // p_k+1 = p_k + omega*(g - B*u_k+1 - C*p_k)
            auto pUpdate = g;
            B.mmv(velocity, pUpdate);
            C.mmv(pressure, pUpdate);
            pUpdate *= w_;
            pressure += pUpdate;
        }
    }

    /*!
     *   \brief Clean up.
     *
     *   \copydoc Preconditioner::post(X&)
     */
    virtual void post (X& x) {}

    //! Category of the preconditioner (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
        return Dune::SolverCategory::sequential;
    }

private:

    void init_(const ParameterTree& params)
    {
        using namespace Dune::Indices;

        if (inexact_)
        {
            linearOperator_ = std::make_shared<LinearOperator>(A_[_0][_0]);
            velocityAMG_ = std::make_unique<VelocityAMG>(linearOperator_, params);

            static const auto determineOmega = Dumux::getParamFromGroup<bool>(paramGroup_, "LinearSolver.PreconditionerDetermineOmega", false);
            if (determineOmega)
                w_ = estimateOmega_();
        }
        else
            velocityUMFPack_ = std::make_unique<VelocityUMFPack>(A_[_0][_0]);
    }

    scalar_field_type estimateOmega_() const
    {
        using namespace Dune::Indices;
        auto& A = A_[_0][_0];
        auto& Bt = A_[_0][_1];
        auto& B = A_[_1][_0];

        VelocityVector x(A.M());
        x = 1.0;

        scalar_field_type omega = 0.0;

        static const auto iterations = Dumux::getParamFromGroup<std::size_t>(paramGroup_, "LinearSolver.PreconditionerPowerLawIterations", 5);

        // apply power iteration x_k+1 = M*x_k/|M*x_k| for the matrix M = -B*Ainv*Bt
        for (std::size_t i = 0; i < iterations; ++i)
        {
            // btx = Bt*x
            VelocityVector btx(x.size());
            Bt.mv(x, btx);

            // ainvbtx = Ainv*(Bt*x)
            auto ainvbtx = x;
            velocityAMG_->pre(ainvbtx, btx);
            velocityAMG_->apply(ainvbtx, btx);
            velocityAMG_->post(ainvbtx);

            // v = M*x = -B*(Ainv*Bt*x)
            VelocityVector v(x.size());
            B.mv(ainvbtx, v);
            v *= -1.0;

            // eigenvalue lambda = xt*M*x/(xt*x) = xt*v/(xt*x);
            auto lambda = x.dot(v)/(x.dot(x));

            // relaxation factor omega = 1/lambda;
            omega = 1.0/lambda;

            // new iterate x = M*x/|M*x| = v/|v|
            x = v;
            x /= v.two_norm();
        }

        std::cout << "relaxation factor " << omega << std::endl;
        return omega;
    }

    //! \brief The matrix we operate on.
    const M& A_;
    //! \brief The number of steps to do in apply
    const int n_;
    //! \brief The relaxation factor to use
    scalar_field_type w_;

    Comm comm_;
    std::shared_ptr<LinearOperator> linearOperator_;
    std::unique_ptr<VelocityAMG> velocityAMG_;
    std::unique_ptr<VelocityUMFPack> velocityUMFPack_;
    const std::string paramGroup_;
    const bool inexact_;
};
DUMUX_REGISTER_PRECONDITIONER("suzawa", Dune::defaultPreconditionerBlockLevelCreator<NewSeqUzawa, 2>());

} // end namespace Dumux

#endif // DUMUX_PRECONDITIONERS_HH
