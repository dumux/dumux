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

template<class M, class X, class Y, int l = 1, int dim = 2>
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
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;
    //! \brief scalar type underlying the field_type
    typedef Dune::Simd::Scalar<field_type> scalar_field_type;

    /*! \brief Constructor.
     *
     *   constructor gets all parameters to operate the prec.
     *   \param mat The matrix to operate on.
     *   \param n The number of iterations to perform.
     *   \param w The relaxation factor.
     */
    NewSeqUzawa (const std::shared_ptr<const AssembledLinearOperator<M,X,Y>>& mat, const ParameterTree& configuration)
    : _A_(mat->getmat()), _n(3), _w(1)
    {
        using namespace Dune::Indices;  // for _0, _1, etc.
        inexact_ = Dumux::getParam<bool>("LinearSolver.InexactVelocitySolver", false);

        if (inexact_)
        {
            Dune::Amg::Parameters params(15, 2000, 1.2, 1.6, Dune::Amg::atOnceAccu);
            params.setDefaultValuesIsotropic(dim);
            params.setDebugLevel(0);

            using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<VelocityMatrix, Dune::Amg::FirstDiagonal>>;
            Criterion criterion(params);

            using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
            SmootherArgs smootherArgs;
            smootherArgs.iterations = 1;
            smootherArgs.relaxationFactor = 1;

            linearOperator_ = std::make_unique<LinearOperator>(_A_[_0][_0]);
            velocityAMG_ = std::make_unique<VelocityAMG>(*linearOperator_, criterion, smootherArgs, comm_);

            static auto determineOmega = Dumux::getParam<bool>("LinearSolver.PreconditionerDetermineOmega", false);
            if (determineOmega)
            {
                auto& A = _A_[_0][_0];
                auto& Bt = _A_[_0][_1];
                auto& B = _A_[_1][_0];

                VelocityVector x(A.M());
                x = 1.0;

                static const auto iterations = Dumux::getParam<std::size_t>("LinearSolver.PreconditionerPowerLawIterations", 5);

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
                    _w = 1.0/lambda;

                    // new iterate x = M*x/|M*x| = v/|v|
                    x = v;
                    x /= v.two_norm();
                }

                std::cout << "relaxation factor " << _w << std::endl;
            }
        }
        else
        {
            velocityUMFPack_ = std::make_unique<VelocityUMFPack>(_A_[_0][_0]);
        }
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

        auto& A = _A_[_0][_0];
        auto& Bt = _A_[_0][_1];
        auto& B = _A_[_1][_0];
        auto& C = _A_[_1][_1];

        const auto& f = d[_0];
        const auto& g = d[_1];
        auto& velocity = v[_0];
        auto& pressure = v[_1];
        for (int i = 0; i < C.N(); ++i)
        {
            if (std::abs(C[i][i].frobenius_norm() - 1.0) < 1e-14)
                pressure[i] = g[i];
        }

        for (int k = 0; k < _n; ++k)
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
            pUpdate *= _w;
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
    //! \brief The matrix we operate on.
    const M& _A_;
    //! \brief The number of steps to do in apply
    int _n;
    //! \brief The relaxation factor to use
    scalar_field_type _w;

    Comm comm_;
    std::unique_ptr<LinearOperator> linearOperator_;
    std::unique_ptr<VelocityAMG> velocityAMG_;
    std::unique_ptr<VelocityUMFPack> velocityUMFPack_;
    bool inexact_;
};
DUMUX_REGISTER_PRECONDITIONER("suzawa", Dune::defaultPreconditionerBlockLevelCreator<NewSeqUzawa, 2>());

} // end namespace Dumux

#endif // DUMUX_PRECONDITIONERS_HH
