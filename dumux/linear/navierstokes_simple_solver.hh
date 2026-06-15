// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Preconditioned iterative solver for the incompressible Navier-Stokes
 *        problem using a SIMPLE block preconditioner (matrix-free Schur complement).
 */
#ifndef DUMUX_LINEAR_NAVIERSTOKES_SIMPLE_SOLVER_HH
#define DUMUX_LINEAR_NAVIERSTOKES_SIMPLE_SOLVER_HH

#include <type_traits>
#include <memory>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include <dune/common/parametertree.hh>
#include <dune/common/indices.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/istl/operators.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>

#if HAVE_MPI
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/novlpschwarz.hh>
#endif

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/common/typetraits/matrix.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/preconditioners.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/parallelmatrixadapter.hh>
#include <dumux/linear/symmetrize_constraints.hh>

// Reuse the proven parallel saddle-point infrastructure (operator, scalar product,
// preconditioner wrapper, representation helpers) from the Stokes solver.
#include <dumux/linear/stokes_solver.hh>

namespace Dumux::Detail {

#if HAVE_MPI
//! Project a vector block to CONSISTENT representation (every copy holds the global value)
//! by summing the additive per-rank partials across ranks. Unlike completeCouplingResult,
//! this does NOT zero the non-owner entries afterwards (leaves owner == copy == global).
template<class Comm, class Block>
void makeConsistentBlock(const Comm& comm, Block& v)
{
    zeroOverlap(comm, v);
    comm.addOwnerCopyToOwnerCopy(v, v);
}
#endif

/*!
 * \ingroup Linear
 * \brief Matrix-free Schur complement operator \f$ \hat{S} = -C D^{-1} B \f$ for SIMPLE.
 *
 * Never assembles \f$\hat{S}\f$. Each apply evaluates the chain of matrix-vector products
 * \f$ p \mapsto -C\,(D^{-1}(B\,p)) \f$, where \f$D\f$ is the block diagonal of the velocity
 * operator \f$A\f$. In parallel, intermediate velocity quantities are completed to a
 * consistent representation across ranks and the pressure result is masked to unique, exactly
 * mirroring the conventions of ParallelStokesLinearOperator so the inner Krylov solve works
 * in the same consistent/unique representation as the outer one.
 *
 * \tparam M the full 2x2 saddle-point matrix type
 * \tparam X the full solution vector type
 */
template<class M, class X>
class SimpleSchurOperator
: public Dune::LinearOperator<std::decay_t<decltype(std::declval<X>()[Dune::Indices::_1])>,
                              std::decay_t<decltype(std::declval<X>()[Dune::Indices::_1])>>
{
    using A = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_0][Dune::Indices::_0])>;
    using ABlock = typename A::block_type;
#if HAVE_MPI
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
#endif
public:
    using U = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_0])>;
    using V = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_1])>;
    using domain_type = V;
    using range_type = V;
    using field_type = typename V::field_type;

    SimpleSchurOperator(const M& matrix, const std::vector<ABlock>& dinv,
                        bool nonOverlapping = false
#if HAVE_MPI
                        , std::shared_ptr<Comm> vComm = nullptr
                        , std::shared_ptr<Comm> pComm = nullptr
#endif
                        )
    : matrix_(matrix), dinv_(dinv), nonOverlapping_(nonOverlapping)
#if HAVE_MPI
    , vComm_(std::move(vComm)), pComm_(std::move(pComm))
#endif
    {}

    //! y = (C D^{-1} B - D11) p  (the SIMPLE Schur complement; positive definite for the
    //! pressure-scaled saddle-point system)
    void apply(const V& p, V& y) const override
    {
        using namespace Dune::Indices;

        V pc(p);
        U bv(matrix_[_0][_0].N());
        bv = 0.0;
#if HAVE_MPI
        if (pComm_)
            pComm_->copyOwnerToAll(pc, pc); // p -> consistent
#endif
        matrix_[_0][_1].mv(pc, bv); // B p (additive in parallel)
#if HAVE_MPI
        if (vComm_)
            makeConsistentBlock(*vComm_, bv); // -> consistent global velocity
#endif
        applyDinv_(bv); // D^{-1} (B p), stays consistent

        y = 0.0;
        matrix_[_1][_0].mv(bv, y);  // C (D^{-1} B p) (additive in parallel)
        matrix_[_1][_1].mmv(pc, y); // y -= D11 p (the (1,1) block; additive in parallel)
#if HAVE_MPI
        if (pComm_)
            completeCouplingResult(*pComm_, y); // complete both terms -> unique pressure result
#endif
    }

    void applyscaleadd(field_type alpha, const V& p, V& y) const override
    {
        V Sp(y); Sp = 0.0;
        apply(p, Sp);
        y.axpy(alpha, Sp);
    }

    Dune::SolverCategory::Category category() const override
    {
#if HAVE_MPI
        if (vComm_)
            return nonOverlapping_ ? Dune::SolverCategory::nonoverlapping
                                   : Dune::SolverCategory::overlapping;
#endif
        return Dune::SolverCategory::sequential;
    }

private:
    //! apply the block-diagonal inverse D^{-1} in place: v_i <- dinv_i * v_i
    void applyDinv_(U& v) const
    {
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            auto tmp = v[i];
            dinv_[i].mv(v[i], tmp);
            v[i] = tmp;
        }
    }

    const M& matrix_;
    const std::vector<ABlock>& dinv_;
    bool nonOverlapping_;
#if HAVE_MPI
    std::shared_ptr<Comm> vComm_, pComm_;
#endif
};

/*!
 * \ingroup Linear
 * \brief Matrix-free Jacobi inner preconditioner for the Schur Krylov solve.
 *
 * Applies the inverse of the (cheaply pre-computed) Schur-complement diagonal:
 * \f$ v_i = \mathrm{diag}(\hat S)_i^{-1}\, d_i \f$ (identity if no diagonal is given).
 * In parallel the result is broadcast owner->copies (copyOwnerToAll) so the correction fed
 * back to the (nonoverlapping/overlapping) Krylov solver is in consistent representation.
 */
template<class V>
class SchurJacobiPreconditioner : public Dune::Preconditioner<V, V>
{
#if HAVE_MPI
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
#endif
public:
    using domain_type = V;
    using range_type = V;
    using field_type = typename V::field_type;

    //! \param invDiag pointer to the inverse Schur diagonal (nullptr -> identity)
    SchurJacobiPreconditioner(const V* invDiag, bool nonOverlapping = false
#if HAVE_MPI
                              , std::shared_ptr<Comm> pComm = nullptr
#endif
                              )
    : invDiag_(invDiag), nonOverlapping_(nonOverlapping)
#if HAVE_MPI
    , pComm_(std::move(pComm))
#endif
    {}

    void pre(V&, V&) override {}

    void apply(V& v, const V& d) override
    {
        if (invDiag_)
            for (std::size_t i = 0; i < v.size(); ++i)
                v[i] = d[i] * (*invDiag_)[i][0];
        else
            v = d;
#if HAVE_MPI
        if (pComm_)
            pComm_->copyOwnerToAll(v, v);
#endif
    }

    void post(V&) override {}

    Dune::SolverCategory::Category category() const override
    {
#if HAVE_MPI
        if (pComm_)
            return nonOverlapping_ ? Dune::SolverCategory::nonoverlapping
                                   : Dune::SolverCategory::overlapping;
#endif
        return Dune::SolverCategory::sequential;
    }

private:
    const V* invDiag_;
    bool nonOverlapping_;
#if HAVE_MPI
    std::shared_ptr<Comm> pComm_;
#endif
};

/*!
 * \ingroup Linear
 * \brief SIMPLE block preconditioner for the Navier-Stokes saddle-point problem
 * \f$
   \begin{pmatrix} A & B \\ C & 0 \end{pmatrix}
   \begin{pmatrix} u \\ p \end{pmatrix} =
   \begin{pmatrix} f \\ g \end{pmatrix}.
 * \f$
 *
 * Given a defect \f$(r_u, r_p)\f$ the preconditioner produces an update \f$(\delta u, \delta p)\f$:
 *   1. velocity predictor: \f$\hat u = A^{-1} r_u\f$ (velocity AMG, one cycle)
 *   2. pressure rhs:       \f$\hat r_p = r_p - C \hat u\f$
 *   3. Schur solve:        \f$\hat S\,\delta p = \hat r_p\f$, \f$\hat S = -C D^{-1} B\f$ (inner GMRES,
 *                          matrix-free SimpleSchurOperator, few iterations)
 *   4. velocity corrector: \f$\delta u = \hat u - D^{-1} B\,\delta p\f$
 *   5. relaxation:         \f$\delta p \mathrel{*}= \omega_p\f$
 *
 * \tparam M the full 2x2 MultiTypeBlockMatrix type
 * \tparam X the domain (update) type
 * \tparam Y the range (defect) type
 */
template<class M, class X, class Y>
class SimplePreconditioner : public Dune::Preconditioner<X, Y>
{
    static_assert(Dumux::isMultiTypeBlockMatrix<M>::value && M::M() == 2 && M::N() == 2,
                  "SimplePreconditioner expects a 2x2 MultiTypeBlockMatrix.");

    using A = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_0][Dune::Indices::_0])>;
    using U = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_0])>;
    using V = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_1])>;
    using ABlock = typename A::block_type;
    static constexpr int dim = ABlock::rows;

    enum class Variant { simple, simplec };

#if HAVE_MPI
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
#endif

public:
    using matrix_type = M;
    using domain_type = X;
    using range_type = Y;
    using field_type = typename X::field_type;
    using scalar_field_type = Dune::Simd::Scalar<field_type>;

    //! Constructor (sequential)
    SimplePreconditioner(
        const std::shared_ptr<const Dune::AssembledLinearOperator<M, X, Y>>& fullLinearOperator,
        const Dune::ParameterTree& params)
    : matrix_(fullLinearOperator->getmat())
    , verbosity_(params.get<int>("verbosity"))
    , paramGroup_(params.get<std::string>("ParameterGroup"))
    {
        init_(params);
    }

#if HAVE_MPI
    //! Constructor (parallel)
    SimplePreconditioner(
        const std::shared_ptr<const Dune::AssembledLinearOperator<M, X, Y>>& fullLinearOperator,
        const Dune::ParameterTree& params,
        std::shared_ptr<Comm> velComm,
        std::shared_ptr<Comm> presComm)
    : matrix_(fullLinearOperator->getmat())
    , verbosity_(params.get<int>("verbosity"))
    , paramGroup_(params.get<std::string>("ParameterGroup"))
    , velComm_(std::move(velComm))
    , presComm_(std::move(presComm))
    {
        init_(params);
    }
#endif

    void pre(X& update, Y& currentDefect) override
    {
        using namespace Dune::Indices;
        preconditionerForA_->pre(update[_0], currentDefect[_0]);
    }

    void apply(X& update, const Y& currentDefect) override
    {
        using namespace Dune::Indices;

        // 1. velocity predictor  uhat = A^{-1} r_u  (one velocity AMG cycle)
        U uhat(matrix_[_0][_0].N()); uhat = 0.0;
        auto ru = currentDefect[_0];
        preconditionerForA_->apply(uhat, ru); // uhat consistent

        // 2. pressure rhs  rhat_p = C uhat - r_p  (matches the positive-definite
        //    Shat = C D^{-1} B - D11 so that Shat dp = rhat_p gives the SIMPLE correction)
        V cu(matrix_[_1][_0].N()); cu = 0.0;
        matrix_[_1][_0].mv(uhat, cu); // C uhat (additive in parallel)
#if HAVE_MPI
        if (presComm_)
            completeCouplingResult(*presComm_, cu); // -> unique
#endif
        V rhatP = cu;
        rhatP -= currentDefect[_1];

        // 3. inner Schur solve  Shat dp = rhat_p
        V dp(matrix_[_1][_1].N()); dp = 0.0;
        solveSchur_(dp, rhatP);

        // 4. velocity corrector  du = uhat - D^{-1} B dp
        V dpc(dp);
#if HAVE_MPI
        if (presComm_)
            presComm_->copyOwnerToAll(dpc, dpc); // dp -> consistent
#endif
        U bdp(matrix_[_0][_0].N()); bdp = 0.0;
        matrix_[_0][_1].mv(dpc, bdp); // B dp (additive)
#if HAVE_MPI
        if (velComm_)
            makeConsistentBlock(*velComm_, bdp); // -> consistent
#endif
        applyDinv_(bdp); // D^{-1} B dp, consistent

        // assemble the update; owner values are correct, the outer parallel wrapper
        // broadcasts owner -> ghost (copyOwnerToAll) so non-owner entries are fixed up.
        update[_0] = uhat;
        update[_0] -= bdp;
        update[_1] = dp;
        update[_1] *= omegaP_;
    }

    void post(X& update) override
    {
        using namespace Dune::Indices;
        preconditionerForA_->post(update[_0]);
    }

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    void init_(const Dune::ParameterTree& params)
    {
        omegaP_ = getParamFromGroup<double>(
            paramGroup_, "LinearSolver.Preconditioner.PressureRelaxation", 0.8);
        // The inner Schur system is solved (near-)exactly so the SIMPLE preconditioner behaves
        // as a fixed linear operator (otherwise a loose/variable inner solve fools the outer
        // Krylov solver). A cheap matrix-free Jacobi (Schur-diagonal) accelerates this.
        schurIterations_ = getParamFromGroup<int>(
            paramGroup_, "LinearSolver.Preconditioner.SchurIterations", 100);
        schurReduction_ = getParamFromGroup<double>(
            paramGroup_, "LinearSolver.Preconditioner.SchurReduction", 1e-6);

        const auto v = getParamFromGroup<std::string>(
            paramGroup_, "LinearSolver.Preconditioner.Variant", "SIMPLE");
        auto lower = v;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        variant_ = (lower == "simplec") ? Variant::simplec : Variant::simple;

        buildDinv_();
        computeSchurDiag_();
        initVelocityPreconditioner_(params);
    }

    //! build the block-diagonal inverse D^{-1} of the velocity operator A
    void buildDinv_()
    {
        using namespace Dune::Indices;
        const auto& A00 = matrix_[_0][_0];
        dinv_.resize(A00.N());
        for (auto row = A00.begin(); row != A00.end(); ++row)
        {
            const auto i = row.index();
            ABlock d(0.0);
            if (variant_ == Variant::simplec)
            {
                // SIMPLEC: D = row sum of A (lumped), captures off-diagonal coupling
                for (auto col = row->begin(); col != row->end(); ++col)
                    d += *col;
            }
            else
            {
                // SIMPLE: D = diagonal block of A
                d = A00[i][i];
            }
            d.invert();
            dinv_[i] = d;
        }
    }

    void applyDinv_(U& v) const
    {
        for (std::size_t i = 0; i < v.size(); ++i)
        {
            auto tmp = v[i];
            dinv_[i].mv(v[i], tmp);
            v[i] = tmp;
        }
    }

    //! cheap matrix-free Schur diagonal via sparse contraction (no matrix-matrix product):
    //! diag(Shat)_i = sum_k C[i][k] D^{-1}_k B[k][i] - D11[i][i]; used as a Jacobi
    //! preconditioner for the inner Schur Krylov solve.
    void computeSchurDiag_()
    {
        using namespace Dune::Indices;
        const auto& C = matrix_[_1][_0];
        const auto& B = matrix_[_0][_1];
        const auto& D11 = matrix_[_1][_1];

        schurInvDiag_.resize(C.N());
        schurInvDiag_ = 0.0;
        for (auto rowC = C.begin(); rowC != C.end(); ++rowC)
        {
            const auto i = rowC.index();
            field_type diag = 0.0;
            for (auto colC = rowC->begin(); colC != rowC->end(); ++colC)
            {
                const auto k = colC.index(); // velocity dof
                const auto bIt = B[k].find(i);
                if (bIt == B[k].end())
                    continue;
                const auto& Cik = *colC;   // FieldMatrix<1,dim>
                const auto& Bki = *bIt;    // FieldMatrix<dim,1>
                Dune::FieldVector<field_type, dim> bcol(0.0), Dt(0.0);
                for (int a = 0; a < dim; ++a) bcol[a] = Bki[a][0];
                dinv_[k].mv(bcol, Dt);     // D^{-1}_k * B[k][i]
                for (int a = 0; a < dim; ++a) diag += Cik[0][a]*Dt[a];
            }
            if (D11.exists(i, i))
                diag -= D11[i][i][0][0];

            using std::abs;
            schurInvDiag_[i] = (abs(diag) > 1e-30) ? 1.0/diag : 1.0;
        }
    }

    //! solve the (matrix-free) Schur system with an inner Jacobi-preconditioned Krylov solve
    void solveSchur_(V& dp, const V& rhs)
    {
        using SchurOp = SimpleSchurOperator<M, X>;
        auto schurOp = std::make_shared<SchurOp>(
            matrix_, dinv_, nonOverlapping_()
#if HAVE_MPI
            , velComm_, presComm_
#endif
        );

        std::shared_ptr<Dune::ScalarProduct<V>> sp;
        auto prec = std::make_shared<SchurJacobiPreconditioner<V>>(
            &schurInvDiag_, nonOverlapping_()
#if HAVE_MPI
            , presComm_
#endif
        );
#if HAVE_MPI
        if (presComm_)
            sp = Dune::createScalarProduct<V>(*presComm_, schurOp->category());
        else
#endif
            sp = std::make_shared<Dune::SeqScalarProduct<V>>();

        Dune::RestartedGMResSolver<V> solver(
            schurOp, sp, prec, schurReduction_, schurIterations_, schurIterations_, /*verbose*/0);

        V b(rhs);
        Dune::InverseOperatorResult res;
        solver.apply(dp, b, res);
    }

    void initVelocityPreconditioner_(const Dune::ParameterTree& params)
    {
        using namespace Dune::Indices;

        if (getParamFromGroup<bool>(paramGroup_, "LinearSolver.DirectSolverForVelocity", false))
        {
#if HAVE_UMFPACK
            directSolverForA_ = std::make_shared<Dune::UMFPack<A>>(matrix_[_0][_0], verbosity_);
            using Wrap = Dune::InverseOperator2Preconditioner<Dune::InverseOperator<U, U>>;
            preconditionerForA_ = std::make_shared<Wrap>(*directSolverForA_);
#else
            DUNE_THROW(Dune::InvalidStateException, "Selected direct solver but UMFPack is not available.");
#endif
        }
        else
        {
#if HAVE_MPI
            if (velComm_)
            {
                const bool nonov = (velComm_->category() == Dune::SolverCategory::nonoverlapping);
                if (nonov)
                {
                    using VelOp = Dune::NonoverlappingSchwarzOperator<A, U, U, Comm>;
                    using Smoother = Dune::NonoverlappingBlockPreconditioner<Comm, Dune::SeqSSOR<A, U, U>>;
                    auto lopV = std::make_shared<VelOp>(matrix_[_0][_0], *velComm_);
                    preconditionerForA_ = std::make_shared<
                        Dune::Amg::AMG<VelOp, U, Smoother, Comm>>(lopV, params, *velComm_);
                }
                else
                {
                    using VelOp = Dune::OverlappingSchwarzOperator<A, U, U, Comm>;
                    using Smoother = Dune::BlockPreconditioner<U, U, Comm, Dune::SeqSSOR<A, U, U>>;
                    auto lopV = std::make_shared<VelOp>(matrix_[_0][_0], *velComm_);
                    preconditionerForA_ = std::make_shared<
                        Dune::Amg::AMG<VelOp, U, Smoother, Comm>>(lopV, params, *velComm_);
                }
            }
            else
#endif
            {
                using VelLinearOperator = Dune::MatrixAdapter<A, U, U>;
                auto lopV = std::make_shared<VelLinearOperator>(matrix_[_0][_0]);
                preconditionerForA_ = std::make_shared<
                    Dune::Amg::AMG<VelLinearOperator, U, Dumux::ParMTSSOR<A, U, U>>>(lopV, params);
            }
        }
    }

    bool nonOverlapping_() const
    {
#if HAVE_MPI
        if (velComm_)
            return velComm_->category() == Dune::SolverCategory::nonoverlapping;
#endif
        return false;
    }

    const M& matrix_;
    const int verbosity_;
    const std::string paramGroup_;

    std::vector<ABlock> dinv_;
    V schurInvDiag_;
    Variant variant_ = Variant::simple;
    double omegaP_ = 0.8;
    int schurIterations_ = 20;
    double schurReduction_ = 1e-2;

    std::shared_ptr<Dune::Preconditioner<U, U>> preconditionerForA_;
    std::shared_ptr<Dune::InverseOperator<U, U>> directSolverForA_;

#if HAVE_MPI
    std::shared_ptr<Comm> velComm_;
    std::shared_ptr<Comm> presComm_;
#endif
};

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Preconditioned iterative solver for the incompressible Navier-Stokes problem
 * \note Uses a SIMPLE block preconditioner (Detail::SimplePreconditioner) with a matrix-free
 *       Schur complement, suited for (instationary) Navier-Stokes at moderate Reynolds numbers.
 * \note Sequential and MPI-parallel (non-overlapping or overlapping decomposition); reuses the
 *       parallel saddle-point operator/scalar-product/wrapper from StokesSolver.
 */
template<class Matrix, class Vector, class VelocityGG, class PressureGG>
class NavierStokesSimpleSolver
: public LinearSolver
{
    using Preconditioner = Detail::SimplePreconditioner<Matrix, Vector, Vector>;

#if HAVE_MPI
    using VTraits = LinearSolverTraits<VelocityGG>;
    using PTraits = LinearSolverTraits<PressureGG>;
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
#endif

public:
    /*!
     * \brief Constructor
     * \param vGridGeometry grid geometry of the velocity discretization
     * \param pGridGeometry grid geometry of the pressure discretization
     * \param dirichletDofs a vector (same size/shape as the rhs) marking Dirichlet dofs with 1.0
     * \param paramGroup group prefix for parameter lookup
     */
    NavierStokesSimpleSolver(std::shared_ptr<const VelocityGG> vGridGeometry,
                             std::shared_ptr<const PressureGG> pGridGeometry,
                             const Vector& dirichletDofs,
                             const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
    , vGridGeometry_(std::move(vGridGeometry))
    , pGridGeometry_(std::move(pGridGeometry))
    , dirichletDofs_(dirichletDofs)
    {
        params_ = LinearSolverParameters<LinearSolverTraits<VelocityGG>>::createParameterTree(this->paramGroup());
        density_ = getParamFromGroup<double>(this->paramGroup(), "Component.LiquidDensity");
        // The SIMPLE preconditioner is a VARIABLE preconditioner (it contains an inner Krylov
        // Schur solve), so a non-flexible Krylov method (plain GMRES, BiCGSTAB) is silently
        // fooled (it converges in the preconditioned norm while the true residual stalls).
        // Flexible GMRES is the only correct outer solver here; coerce gmres -> fgmres.
        solverType_ = getParamFromGroup<std::string>(this->paramGroup(), "LinearSolver.Type", "fgmres");
        if (solverType_ == "gmres")
            solverType_ = "fgmres";

#if HAVE_MPI
        if (vGridGeometry_->gridView().comm().size() > 1)
            initParallelInfrastructure_();
        else
            scalarProduct_ = std::make_shared<Dune::ScalarProduct<Vector>>();
#else
        scalarProduct_ = std::make_shared<Dune::ScalarProduct<Vector>>();
#endif
    }

    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        auto bTmp = b;
        auto ATmp = A;
        return applyIterativeSolver_(ATmp, x, bTmp);
    }

    Scalar norm(const Vector& b) const
    {
#if HAVE_MPI
        using namespace Dune::Indices;
        if (vComm_ && pComm_ && vGridGeometry_->gridView().comm().size() > 1)
        {
            Vector bc(b);
            vComm_->addOwnerCopyToOwnerCopy(bc[_0], bc[_0]);
            pComm_->addOwnerCopyToOwnerCopy(bc[_1], bc[_1]);
            return scalarProduct_->norm(bc);
        }
#endif
        return scalarProduct_->norm(b);
    }

    std::string name() const
    {
        return "SIMPLE-preconditioned Navier-Stokes solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
#if HAVE_MPI
    void initParallelInfrastructure_()
    {
        isNonOverlapping_ = VTraits::isNonOverlapping(vGridGeometry_->gridView())
                         && PTraits::isNonOverlapping(pGridGeometry_->gridView());

        const auto cat = isNonOverlapping_ ? Dune::SolverCategory::nonoverlapping
                                           : Dune::SolverCategory::overlapping;

        vParallelHelper_ = std::make_shared<ParallelISTLHelper<VTraits>>(
            vGridGeometry_->gridView(), VTraits::dofMapper(*vGridGeometry_));
        pParallelHelper_ = std::make_shared<ParallelISTLHelper<PTraits>>(
            pGridGeometry_->gridView(), PTraits::dofMapper(*pGridGeometry_));

        vComm_ = std::make_shared<Comm>(vGridGeometry_->gridView().comm(), cat);
        vParallelHelper_->createParallelIndexSet(*vComm_);

        pComm_ = std::make_shared<Comm>(pGridGeometry_->gridView().comm(), cat);
        pParallelHelper_->createParallelIndexSet(*pComm_);

        scalarProduct_ = std::make_shared<Detail::ParallelStokesScalarProduct<Vector>>(
            vComm_, pComm_, isNonOverlapping_);
    }
#endif

    bool applyIterativeSolver_(Matrix& A, Vector& x, Vector& b)
    {
        // make Dirichlet boundary conditions symmetric
        if (getParamFromGroup<bool>(this->paramGroup(), "LinearSolver.SymmetrizeDirichlet", true))
            symmetrizeConstraints(A, b, dirichletDofs_);

        // scale the pressure row (block-scale symmetry), same as the Stokes solver
        using namespace Dune::Indices;
        A[_1] *= -1.0/density_;
        b[_1] *= -1.0/density_;

#if HAVE_MPI
        if (vGridGeometry_->gridView().comm().size() > 1)
            return solveParallel_(A, x, b);
#endif
        return solveSequential_(A, x, b);
    }

    bool solveSequential_(Matrix& A, Vector& x, Vector& b)
    {
        auto op = std::make_shared<Dumux::ParallelMultiTypeMatrixAdapter<Matrix, Vector, Vector>>(A);
        auto preconditioner = std::make_shared<Preconditioner>(op, params_.sub("preconditioner"));
        return runSolver_(op, scalarProduct_, preconditioner, x, b);
    }

#if HAVE_MPI
    bool solveParallel_(Matrix& A, Vector& x, Vector& b)
    {
        using namespace Dune::Indices;

        if (isNonOverlapping_)
            prepareParallelLinearSystem_(A, b);

        auto innerOp = std::make_shared<Dumux::ParallelMultiTypeMatrixAdapter<Matrix, Vector, Vector>>(A);
        auto op = std::make_shared<Detail::ParallelStokesLinearOperator<Matrix, Vector, Vector>>(
            A, vComm_, pComm_, isNonOverlapping_);
        auto seqPrec = std::make_shared<Preconditioner>(innerOp, params_.sub("preconditioner"), vComm_, pComm_);
        auto prec = std::make_shared<Detail::ParallelStokesPreconditioner<Preconditioner, Vector, Vector>>(
            seqPrec, vComm_, pComm_, isNonOverlapping_);

        const bool converged = runSolver_(op, scalarProduct_, prec, x, b);

        vComm_->copyOwnerToAll(x[_0], x[_0]);
        pComm_->copyOwnerToAll(x[_1], x[_1]);

        return converged;
    }

    void prepareParallelLinearSystem_(Matrix& A, Vector& b)
    {
        using namespace Dune::Indices;
        using GridView = typename VelocityGG::GridView;
        static constexpr int dim = GridView::dimension;
        static constexpr std::size_t numCodims = dim + 1;

        const auto& gv = vGridGeometry_->gridView();
        const auto& velMapper = VTraits::dofMapper(*vGridGeometry_);
        const auto velCodims = activeCodimsBitset_<VTraits, numCodims>();

        // Velocity diagonal block A00: extendMatrix + sumEntries (summed form needed by the
        // velocity AMG / NonoverlappingSchwarzOperator and by D^{-1}). Off-diagonals B, C are
        // kept RAW (additive) and completed at the mat-vec result level.
        using VelBlock = std::decay_t<decltype(A[_0][_0])>;
        MultiCodimParallelMatrixHelper<VelBlock, GridView, std::decay_t<decltype(velMapper)>, numCodims>
            velHelper(gv, velMapper, velCodims);
        velHelper.extendMatrix(A[_0][_0], [](auto){ return false; });
        velHelper.sumEntries(A[_0][_0]);

        // Make the rhs unique (additive assembly summed across ranks, then mask to unique).
        using VelVec = std::decay_t<decltype(b[_0])>;
        using PresVec = std::decay_t<decltype(b[_1])>;
        using VNonoverlapping = typename VTraits::template ParallelNonoverlapping<VelBlock, VelVec>;
        using PNonoverlapping = typename PTraits::template ParallelNonoverlapping<
            std::decay_t<decltype(A[_1][_1])>, PresVec>;
        prepareVectorParallel<VTraits, VNonoverlapping>(b[_0], *vParallelHelper_);
        prepareVectorParallel<PTraits, PNonoverlapping>(b[_1], *pParallelHelper_);
        Detail::makeUnique(*vComm_, b[_0]);
        Detail::makeUnique(*pComm_, b[_1]);
    }

    template<class Traits, std::size_t numCodims>
    static std::bitset<numCodims> activeCodimsBitset_()
    {
        if constexpr (requires { Traits::dofCodims; })
            return Traits::dofCodims;
        else
        {
            std::bitset<numCodims> bits;
            bits.set(Traits::dofCodim);
            return bits;
        }
    }
#endif // HAVE_MPI

    template<class LinearOperator, class SP, class Prec>
    bool runSolver_(std::shared_ptr<LinearOperator> op,
                    std::shared_ptr<SP> sp,
                    std::shared_ptr<Prec> prec,
                    Vector& x, Vector& b)
    {
        // The SIMPLE preconditioner contains an inner Krylov (Schur) solve, which makes it a
        // VARIABLE (non-constant) preconditioner. A standard GMRES assumes a fixed linear
        // preconditioner and is fooled by this (its residual recurrence diverges from the true
        // residual). A flexible solver (FGMRES) stores the preconditioned vectors explicitly and
        // converges on the true residual, so it is the correct/default choice here.
        std::unique_ptr<Dune::InverseOperator<Vector, Vector>> solver;
        if (solverType_ == "bicgstab")
            solver = std::make_unique<Dune::BiCGSTABSolver<Vector>>(op, sp, prec, params_);
        else if (solverType_ == "gmres")
            solver = std::make_unique<Dune::RestartedGMResSolver<Vector>>(op, sp, prec, params_);
        else if (solverType_ == "fgmres")
            solver = std::make_unique<Dune::RestartedFlexibleGMResSolver<Vector>>(op, sp, prec, params_);
        else
            DUNE_THROW(Dune::NotImplemented, "Solver choice " << solverType_ << " is not implemented");

        solver->apply(x, b, result_);
        return result_.converged;
    }

    double density_;
    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
    std::shared_ptr<const VelocityGG> vGridGeometry_;
    std::shared_ptr<const PressureGG> pGridGeometry_;
    const Vector& dirichletDofs_;
    std::string solverType_;

    std::shared_ptr<Dune::ScalarProduct<Vector>> scalarProduct_;

#if HAVE_MPI
    std::shared_ptr<ParallelISTLHelper<VTraits>> vParallelHelper_;
    std::shared_ptr<ParallelISTLHelper<PTraits>> pParallelHelper_;
    std::shared_ptr<Comm> vComm_, pComm_;
    bool isNonOverlapping_ = true;
#endif
};

} // end namespace Dumux

#endif
