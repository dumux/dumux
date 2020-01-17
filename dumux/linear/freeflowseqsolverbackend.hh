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
 * \brief Dumux sequential linear solver backends
 */
#ifndef DUMUX_FREEFLOW_SEQ_SOLVER_BACKEND_HH
#define DUMUX_FREEFLOW_SEQ_SOLVER_BACKEND_HH

#include <type_traits>
#include <tuple>
#include <utility>

#include "seqsolverbackend.hh"
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/umfpack.hh>
#include <dune/common/version.hh>
#include <dune/common/hybridutilities.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/matrix.hh>
#include <dumux/common/typetraits/utility.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/amgbackend.hh>


namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A simple ilu0 block diagonal preconditioner
 */
template<class M, class X, class Y, template<class Mat, class Vec, size_t i> class Smoother, int blockLevel = 2>
class FreeFlowBlockDiagAMGPreconditioner : public Dune::Preconditioner<X, Y>
{
    template<std::size_t i>
    using DiagBlockType = std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;

    // template<std::size_t i, std::size_t j>
    // using OffDiagBlockType = std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<j>{}])>;

    template<std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<X>()[Dune::index_constant<i>{}])>;

    // template<std::size_t i>
    // using Smoother = Dune::SeqSSOR<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>>;

    template<std::size_t i>
    using LinearOperator = Dune::MatrixAdapter<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>>;

    template<std::size_t i>
    using ScalarProduct = Dune::SeqScalarProduct<VecBlockType<i>>;

    template<std::size_t i>
    using BlockAMG = Dune::Amg::AMG<LinearOperator<i>, VecBlockType<i>, Smoother<M, X, i>>;

    // template<std::size_t i>
    // using BlockILU = Dune::SeqILU<DiagBlockType<i>, VecBlockType<i>, VecBlockType<i>, blockLevel-1>;

    template<std::size_t i>
    using BlockIdentity = Dune::Richardson<VecBlockType<i>, VecBlockType<i>>;

    template<std::size_t i>
    using BlockZero = Dune::Richardson<VecBlockType<i>,VecBlockType<i>>;

    // the actual types of preconditioners
    using VelocityBlockAMG = BlockAMG<1>;
    using PressureBlockIdentity = BlockIdentity<0>;
    using VelocityPressureCouplingBlockPreconditioner = BlockZero<1>;
    using PressureVelocityCouplingBlockPreconditioner = BlockZero<0>;

    using FirstRowPrecMat = Dune::MultiTypeBlockVector<PressureBlockIdentity, PressureVelocityCouplingBlockPreconditioner>;
    using SecondRowPrecMat = Dune::MultiTypeBlockVector<VelocityPressureCouplingBlockPreconditioner, VelocityBlockAMG>;
    using PreconditionerMatrix = Dune::MultiTypeBlockVector<FirstRowPrecMat, SecondRowPrecMat>;

public:
    //! \brief The matrix type the preconditioner is for.
    using matrix_type = typename std::decay_t<M>;
    //! \brief The domain type of the preconditioner.
    using domain_type = X;
    //! \brief The range type of the preconditioner.
    using range_type = Y;
    //! \brief The field type of the preconditioner.
    using field_type = typename X::field_type;

    // using Dune::Indices;

    /*! \brief Constructor.

       Constructor gets all parameters to operate the prec.
       \param lop The linear operator
       \param c The criterion
       \param sa The smoother arguments
     */
    template<class LOP, class Criterion, class SmootherArgs>
    FreeFlowBlockDiagAMGPreconditioner(const LOP& lop, const Criterion& c, const SmootherArgs& sa, const M& m, int precondIter, double relaxation)
    : preconditionerMatrix_(FirstRowPrecMat(PressureBlockIdentity(1.0), PressureVelocityCouplingBlockPreconditioner(0.0)),
                            SecondRowPrecMat(VelocityPressureCouplingBlockPreconditioner(0.0), VelocityBlockAMG(lop, c, sa)))
    {
        static_assert(blockLevel >= 2, "Only makes sense for MultiTypeBlockMatrix!");
    }

    /*!
     * \brief Prepare the preconditioner.
     *
     * A solver solves a linear operator equation A(v)=d by applying
     * one or several steps of the preconditioner. The method pre()
     * is called before the first apply operation.
     * d and v are right hand side and solution vector of the linear
     * system respectively. It may. e.g., scale the system, allocate memory or
     * compute a (I)LU decomposition.
     * Note: The ILU decomposition could also be computed in the constructor
     * or with a separate method of the derived method if several
     * linear systems with the same matrix are to be solved.
     *
     * \param v The left hand side of the equation.
     * \param d The right hand side of the equation.
     */
    void pre (X& v, Y& d) final
    {
        // using namespace Dune::Indices;
        // preconditionerMatrix_[_1][_1].pre(v[_1], d[_1]);


        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(preconditionerMatrix_)), [&](const auto i)
        {
            forEach(integralRange(Dune::Hybrid::size(preconditionerMatrix_)), [&](const auto j)
            {
                if (i == 1) // TODO ???
                    preconditionerMatrix_[i][j].pre(v[i], d[i]);
            });
        });
    }

    /*!
     * \brief Apply one step of the preconditioner to the system A(v)=d.
     *
     * On entry v=0 and d=b-A(x) (although this might not be
     * computed in that way. On exit v contains the update, i.e
     * one step computes \f$ v = M^{-1} d \f$ where \f$ M \f$ is the
     * approximate inverse of the operator \f$ A \f$ characterizing
     * the preconditioner.
     * \param v The update to be computed
     * \param d The current defect.
     */
    void apply (X& v, const Y& d) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(preconditionerMatrix_)), [&](const auto i)
        {
            forEach(integralRange(Dune::Hybrid::size(preconditionerMatrix_)), [&](const auto j)
            {
                if (i == j) // TODO why does applying the off-diagonal precondioners not work (i.e. lead to wrong results)?
                preconditionerMatrix_[i][j].apply(v[i], d[i]);
            });
        });

        // // does the same as
        // using namespace Dune::Indices;
        // preconditionerMatrix_[_1][_1].apply(v[_1], d[_1]);
        // preconditionerMatrix_[_0][_0].apply(v[_0], d[_0]); // equals v[_0] = d[_0]
        // preconditionerMatrix_[_0][_1].apply(v[_0], d[_0]);
        // preconditionerMatrix_[_1][_0].apply(v[_1], d[_1]);
        // v[_0] = d[_0];
    }

    /*!
     * \brief Clean up.
     *
     * This method is called after the last apply call for the
     * linear system to be solved. Memory may be deallocated safely
     * here. v is the solution of the linear equation.
     *
     * \param v The right hand side of the equation.
     */
    void post (X& v) final
    {
        using namespace Dune::Hybrid;
        forEach(integralRange(Dune::Hybrid::size(preconditionerMatrix_)), [&](const auto i)
        {
            forEach(integralRange(Dune::Hybrid::size(preconditionerMatrix_)), [&](const auto j)
            {
                if (i == j) // TODO: ???
                preconditionerMatrix_[i][j].post(v[i]);
            });
        });
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const final
    {
        return Dune::SolverCategory::sequential;
    }

private:

    PreconditionerMatrix preconditionerMatrix_;
};

/*!
 * \ingroup Linear
 * \brief A simple ilu0 block diagonal preconditioned BiCGSTABSolver
 * \note expects a system as a multi-type block-matrix
 * | A  B |
 * | C  D |
 */
class BlockDiagAMGRestartedGMRESSolver : public LinearSolver
{
    template<class M, std::size_t i>
    using DiagBlockType = std::decay_t<decltype(std::declval<M>()[Dune::index_constant<i>{}][Dune::index_constant<i>{}])>;

    template<class X, std::size_t i>
    using VecBlockType = std::decay_t<decltype(std::declval<X>()[Dune::index_constant<i>{}])>;

    template<class M, class X, std::size_t i>
    using Smoother = Dune::SeqSSOR<DiagBlockType<M, i>, VecBlockType<X, i>, VecBlockType<X, i>>;
    // using Smoother = Dune::SeqSOR<DiagBlockType<M, i>, VecBlockType<X, i>, VecBlockType<X, i>>;
    // using Smoother = Dune::SeqILU<DiagBlockType<M, i>, VecBlockType<X, i>, VecBlockType<X, i>, /*blockLevel*/2-1>;

    template<class M, class X, std::size_t i>
    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother<M, X, i>>::Arguments;

    template<class M, std::size_t i>
    using Criterion = Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<DiagBlockType<M, i>, Dune::Amg::FirstDiagonal>>;

    template<class M, class X, std::size_t i>
    using LinearOperator = Dune::MatrixAdapter<DiagBlockType<M, i>, VecBlockType<X, i>, VecBlockType<X, i>>;

public:
    using LinearSolver::LinearSolver;

    // Solve saddle-point problem using a Schur complement based preconditioner
    template<int precondBlockLevel = 2, class Matrix, class Vector>
    bool solve(const Matrix& m, Vector& x, const Vector& b)
    {
        //! \todo Check whether the default accumulation mode atOnceAccu is needed.
        static const int dim = getParam<int>("LinearSolver.AMG.Dimension", 2);
        static const int maxLevel = getParam<int>("LinearSolver.AMG.MaxLevel", 15);
        static const int coarsenTarget = getParam<int>("LinearSolver.AMG.CoarsenTarget", 2000);
        static const double minCoarsenRate = getParam<double>("LinearSolver.AMG.MinCoarsenRate", 1.6);
        static const double prolongDamp = getParam<double>("LinearSolver.AMG.ProlongDamp", 1.2);
        static const int gamma = getParam<int>("LinearSolver.AMG.Gamma", 1);
        static const int preSmoothSteps = getParam<int>("LinearSolver.AMG.PreSmoothSteps", 2);
        static const int postSmoothSteps = getParam<int>("LinearSolver.AMG.PostSmoothSteps", 2);
        Dune::Amg::Parameters params(maxLevel, coarsenTarget, minCoarsenRate, prolongDamp, Dune::Amg::atOnceAccu);
        params.setDefaultValuesIsotropic(dim);
        params.setDebugLevel(this->verbosity());
        params.setGamma(gamma);
        params.setNoPreSmoothSteps(preSmoothSteps);
        params.setNoPostSmoothSteps(postSmoothSteps);


        auto criterion = Criterion<Matrix, 1>(params);
        auto smootherArgs = SmootherArgs<Matrix, Vector, 1>();
        smootherArgs.iterations = this->precondIter();
        smootherArgs.relaxationFactor = this->relaxation();
        auto linearOperator = LinearOperator<Matrix, Vector, 1>(m[Dune::index_constant<1>{}][Dune::index_constant<1>{}]);
        FreeFlowBlockDiagAMGPreconditioner<Matrix, Vector, Vector, Smoother> preconditioner(linearOperator, criterion, smootherArgs, m, this->precondIter(), this->relaxation());
        Dune::MatrixAdapter<Matrix, Vector, Vector> op(m);

        auto bTmp = b;

        // calculate reduction criterion
        op.apply(x, bTmp);
        bTmp -= b;
        const auto resNorm = bTmp.two_norm();
        const auto reduction = std::max(1e-16 / resNorm, this->residReduction());

        bTmp = b;

        const int restartGMRes = getParam<int>("LinearSolver.GMResRestart");
        Dune::RestartedGMResSolver<Vector> solver(op, preconditioner, /*this->residReduction()*/reduction, restartGMRes,
                                                  this->maxIter(), this->verbosity());
            // Dune::BiCGSTABSolver<Vector> solver(op, preconditioner, /*this->residReduction()*/reduction,
            //                                           this->maxIter(), this->verbosity());
            // Dune::GeneralizedPCGSolver<Vector> solver(op, preconditioner, /*this->residReduction()*/reduction,
            //                                           this->maxIter(), this->verbosity());
        solver.apply(x, bTmp, result_);

        return result_.converged;
    }

    const Dune::InverseOperatorResult& result() const
    {
      return result_;
    }

    std::string name() const
    { return "block-diagonal AMG preconditioned RestartedGMResSolver solver"; }

private:

    Dune::InverseOperatorResult result_;
};

} // end namespace Dumux

#endif
