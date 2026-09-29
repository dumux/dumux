// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Preconditioned iterative solver for the incompressible Stokes problem
 */
#ifndef DUMUX_LINEAR_STOKES_SOLVER_HH
#define DUMUX_LINEAR_STOKES_SOLVER_HH

#include <type_traits>
#include <memory>
#include <tuple>
#include <map>
#include <algorithm>
#include <cmath>

#include <dune/common/parametertree.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/exceptions.hh>

#include <dune/istl/matrixindexset.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/paamg/amg.hh>

#if HAVE_MPI
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/novlpschwarz.hh>
#endif

#include <dumux/common/math.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/preconditioners.hh>
#include <dumux/linear/linearsolverparameters.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/parallelmatrixadapter.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/linear/symmetrize_constraints.hh>
#include <dumux/assembly/jacobianpattern.hh>

namespace Dumux::Detail {

#if HAVE_MPI
//! Keep the entries of owned dofs and set those of all other shared dofs to zero,
//! i.e. the unique representation (Blatt & Bastian §2.5)
template<class Comm, class Block>
void makeUnique(const Comm& comm, Block& v)
{
    using Attr = Dune::OwnerOverlapCopyAttributeSet;
    for (const auto& pair : comm.indexSet())
        if (pair.local().attribute() != Attr::owner)
            v[pair.local().local()] = 0.0;
}

//! Set the entries of overlap dofs to zero
template<class Comm, class Block>
void zeroOverlap(const Comm& comm, Block& v)
{
    using Attr = Dune::OwnerOverlapCopyAttributeSet;
    for (const auto& pair : comm.indexSet())
        if (pair.local().attribute() == Attr::overlap)
            v[pair.local().local()] = 0.0;
}

/*!
 * \brief Complete the product y = M x of an off-diagonal block M with a consistent x to the unique representation
 * \note On a non-overlapping decomposition M is stored additively. Summing the partial products of all processes
 *       also completes the rows at the processor border whose stencil reaches dofs that only a neighbouring
 *       process assembles.
 */
template<class Comm, class Block>
void completeCouplingResult(const Comm& comm, Block& y)
{
    // on an overlapping decomposition the owned rows are complete already
    if (comm.category() == Dune::SolverCategory::overlapping)
    {
        makeUnique(comm, y);
        return;
    }

    zeroOverlap(comm, y);
    comm.addOwnerCopyToOwnerCopy(y, y);
    makeUnique(comm, y);
}
#endif

/*!
 * \ingroup Linear
 * \brief A Stokes preconditioner (saddle-point problem) for the problem
 * \f$
 \begin{pmatrix} A & B \\ C & 0 \end{pmatrix}
 \begin{pmatrix} u \\ p \end{pmatrix} =
 \begin{pmatrix} f \\ g \end{pmatrix},
 * \f$
 *
 * where A is the discrete velocity operator and B and C are discrete gradient and divergence operators.
 * This preconditioner is especially suited for solving the incompressible Stokes equations.
 *
 * \tparam M Type of the matrix.
 * \tparam X Type of the update.
 * \tparam Y Type of the defect.
 * \tparam l Preconditioner block level
 */
template<class M, class X, class Y, int l = 2>
class StokesPreconditioner : public Dune::Preconditioner<X,Y>
{
    static_assert(Dumux::isMultiTypeBlockMatrix<M>::value && M::M() == 2 && M::N() == 2, "Expects a 2x2 MultiTypeBlockMatrix.");
    static_assert(l == 2, "StokesPreconditioner expects a block level of 2 for Matrix type M.");

    using A = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_0][Dune::Indices::_0])>;
    using U = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_0])>;

    using P = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_1][Dune::Indices::_1])>;
    using V = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_1])>;

    enum class Mode { symmetric, triangular,  diagonal };

#if HAVE_MPI
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
#endif

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
    //! \brief the type of the pressure operator
    using PressureLinearOperator = Dune::MatrixAdapter<P,V,V>;

    /*!
     * \brief Constructor (sequential)
     * \param fullLinearOperator the Stokes linear operator
     * \param pressureLinearOperator the linear operator for the pressure space preconditioner
     * \param params a parameter tree for the preconditioner configuration
     */
    StokesPreconditioner(
        const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& fullLinearOperator,
        const std::shared_ptr<const PressureLinearOperator>& pressureLinearOperator,
        const Dune::ParameterTree& params
    )
    : matrix_(fullLinearOperator->getmat())
    , pmatrix_(pressureLinearOperator->getmat())
    , verbosity_(params.get<int>("verbosity"))
    , paramGroup_(params.get<std::string>("ParameterGroup"))
    {
        initMode_(params);
        initPreconditioner_(params);
    }

#if HAVE_MPI
    /*!
     * \brief Constructor (parallel)
     * \param fullLinearOperator the Stokes linear operator
     * \param pressureLinearOperator the linear operator for the pressure space preconditioner
     * \param params a parameter tree for the preconditioner configuration
     * \param velComm OwnerOverlapCopy communication for the velocity DOFs
     * \param presComm OwnerOverlapCopy communication for the pressure DOFs
     */
    StokesPreconditioner(
        const std::shared_ptr<const Dune::AssembledLinearOperator<M,X,Y>>& fullLinearOperator,
        const std::shared_ptr<const PressureLinearOperator>& pressureLinearOperator,
        const Dune::ParameterTree& params,
        std::shared_ptr<Comm> velComm,
        std::shared_ptr<Comm> presComm = nullptr
    )
    : matrix_(fullLinearOperator->getmat())
    , pmatrix_(pressureLinearOperator->getmat())
    , verbosity_(params.get<int>("verbosity"))
    , paramGroup_(params.get<std::string>("ParameterGroup"))
    , velComm_(std::move(velComm))
    , presComm_(std::move(presComm))
    {
        initMode_(params);
        initPreconditioner_(params);
    }
#endif

    /*!
     * \brief Prepare the preconditioners of the velocity and the pressure block
     */
    void pre(X& update, Y& currentDefect) override
    {
        using namespace Dune::Indices;
        preconditionerForA_->pre(update[_0], currentDefect[_0]);
        preconditionerForP_->pre(update[_1], currentDefect[_1]);
    }

    /*!
     * \brief Apply the preconditioner
     *
     * \param update The update to be computed.
     * \param currentDefect The current defect.
     *
     * The currentDefect has be be in a consistent representation,
     * Definition 2.3 Blatt and Bastian (2009) https://doi.org/10.1504/IJCSE.2008.021112
     * The update is initially zero. At exit the update has to be
     * in a consistent representation. This usually requires communication.
     */
    void apply(X& update, const Y& currentDefect) override
    {
        using namespace Dune::Indices;

        if (mode_ == Mode::symmetric)
        {
            update = applyRightBlock_(currentDefect);
            update = applyDiagBlock_(update);
            update = applyLeftBlock_(update);
        }
        else if (mode_ == Mode::triangular)
        {
            update = applyRightBlock_(currentDefect);
            update = applyDiagBlock_(update);
        }
        else
            update = applyDiagBlock_(currentDefect);
    }

    /*!
     * \brief Clean up the preconditioners of the velocity and the pressure block
     */
    void post(X& update) override
    {
        using namespace Dune::Indices;
        preconditionerForA_->post(update[_0]);
        preconditionerForP_->post(update[_1]);
    }

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    void initMode_(const Dune::ParameterTree& params)
    {
        const auto mode = getParamFromGroup<std::string>(
            paramGroup_, "LinearSolver.Preconditioner.Mode", "Diagonal"
        );

        auto lower = mode;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower == "symmetric")
            mode_ = Mode::symmetric;
        else if (lower == "triangular")
            mode_ = Mode::triangular;
        else
            mode_ = Mode::diagonal;
    }

    X applyRightBlock_(const Y& d)
    {
        using namespace Dune::Indices;

        // right bit of LDU decomposition
        // applied to d
        //
        // | I         0 |
        // | -C*inv(A) I |
        //

        auto dTmp0 = d[_0];
        auto vTmp = d; vTmp = 0.0;

        // invert velocity block (apply to first block of d), the result is consistent
        applyPreconditionerForA_(vTmp[_0], dTmp0);

        // then multiply with C, in parallel completed to the unique representation of d[_1]
        matrix_[_1][_0].mv(vTmp[_0], vTmp[_1]);
#if HAVE_MPI
        if (presComm_)
            completeCouplingResult(*presComm_, vTmp[_1]);
#endif

        // and subtract from d
        auto v = d;
        v[_0] = vTmp[_0]; // already do A^-1 d of the diagonal block because we already computed it here
        v[_1] -= vTmp[_1];
        return v;
    }

    X applyDiagBlock_(const Y& d)
    {
        using namespace Dune::Indices;

        // diagonal middle bit of LDU decomposition
        auto dTmp = d;
        auto v = d; v = 0.0;

        // invert velocity block
        if (mode_ == Mode::diagonal)
            applyPreconditionerForA_(v[_0], dTmp[_0]);

        // reuse the already computed A^-1 d (see applyRightBlock_)
        else
            v[_0] = dTmp[_0];

        // invert pressure block
        applyPreconditionerForP_(v[_1], dTmp[_1]);
        if (pressureDiagonal_)
            for (std::size_t i = 0; i < v[_1].size(); ++i)
                v[_1][i] += (*pressureDiagonal_)[i]*d[_1][i];

        return v;
    }

    X applyLeftBlock_(const Y& d)
    {
        using namespace Dune::Indices;

        // left bit of LDU decomposition
        // applied to d
        //
        // | I  -inv(A)*B |
        // | 0       I    |
        //

        auto dTmp = d;
        auto vTmp = d; vTmp = 0.0;

        // multiply with B, in parallel completed to the unique defect the velocity preconditioner expects
        matrix_[_0][_1].umv(dTmp[_1], vTmp[_0]);
#if HAVE_MPI
        if (velComm_)
            completeCouplingResult(*velComm_, vTmp[_0]);
#endif

        // invert velocity block (apply to first block of d)
        auto vTmp0 = vTmp[_0]; vTmp0 = 0.0;
        applyPreconditionerForA_(vTmp0, vTmp[_0]);

        // and subtract from d
        auto v = d;
        v[_0] -= vTmp0;

        return v;
    }

    void initPreconditioner_(const Dune::ParameterTree& params)
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
        // a few smoothing sweeps suffice where the velocity block is dominated by its storage term
        else if (getParamFromGroup<bool>(paramGroup_, "LinearSolver.Preconditioner.SSORForVelocity", false))
        {
            const auto sweeps = getParamFromGroup<int>(paramGroup_, "LinearSolver.Preconditioner.VelocitySweeps", 2);
            const auto relaxation = getParamFromGroup<scalar_field_type>(paramGroup_, "LinearSolver.Preconditioner.VelocityRelaxation", 1.0);
            auto ssor = std::make_shared<Dune::SeqSSOR<A,U,U>>(matrix_[_0][_0], sweeps, relaxation);
#if HAVE_MPI
            if (velComm_)
            {
                if (velComm_->category() == Dune::SolverCategory::nonoverlapping)
                    preconditionerForA_ = std::make_shared<Dune::NonoverlappingBlockPreconditioner<Comm, Dune::SeqSSOR<A,U,U>>>(ssor, *velComm_);
                else
                    preconditionerForA_ = std::make_shared<Dune::BlockPreconditioner<U, U, Comm, Dune::SeqSSOR<A,U,U>>>(ssor, *velComm_);
            }
            else
#endif
                preconditionerForA_ = ssor;
        }
        else
        {
#if HAVE_MPI
            if (velComm_)
            {
                const bool nonov = (velComm_->category()
                                    == Dune::SolverCategory::nonoverlapping);
                if (nonov)
                {
                    using VelOp = Dune::NonoverlappingSchwarzOperator<A, U, U, Comm>;
                    using Smoother = Dune::NonoverlappingBlockPreconditioner<Comm, Dune::SeqSSOR<A,U,U>>;
                    auto lopV = std::make_shared<VelOp>(matrix_[_0][_0], *velComm_);
                    preconditionerForA_ = std::make_shared<
                        Dune::Amg::AMG<VelOp, U, Smoother, Comm>>(lopV, params, *velComm_);
                }
                else
                {
                    using VelOp = Dune::OverlappingSchwarzOperator<A, U, U, Comm>;
                    using Smoother = Dune::BlockPreconditioner<U, U, Comm, Dune::SeqSSOR<A,U,U>>;
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
                    Dune::Amg::AMG<VelLinearOperator, U, Dumux::ParMTSSOR<A,U,U>>
                >(lopV, params);
            }
        }

        if (getParamFromGroup<bool>(paramGroup_, "LinearSolver.DirectSolverForPressure", false))
        {
#if HAVE_UMFPACK
            directSolverForP_ = std::make_shared<Dune::UMFPack<P>>(pmatrix_, verbosity_);
            using Wrap = Dune::InverseOperator2Preconditioner<Dune::InverseOperator<V, V>>;
            preconditionerForP_ = std::make_shared<Wrap>(*directSolverForP_);
#else
            DUNE_THROW(Dune::InvalidStateException, "Selected direct solver but UMFPack is not available.");
#endif
        }
        // a Poisson-type pressure operator, as for small time steps, needs more than a Jacobi iteration
        else if (getParamFromGroup<bool>(paramGroup_, "LinearSolver.Preconditioner.AMGForPressure", false))
        {
#if HAVE_MPI
            if (presComm_)
            {
                if (presComm_->category() == Dune::SolverCategory::nonoverlapping)
                {
                    using PressOp = Dune::NonoverlappingSchwarzOperator<P, V, V, Comm>;
                    using Smoother = Dune::NonoverlappingBlockPreconditioner<Comm, Dune::SeqSSOR<P,V,V>>;
                    auto lopP = std::make_shared<PressOp>(pmatrix_, *presComm_);
                    preconditionerForP_ = std::make_shared<
                        Dune::Amg::AMG<PressOp, V, Smoother, Comm>>(lopP, params, *presComm_);
                }
                else
                {
                    using PressOp = Dune::OverlappingSchwarzOperator<P, V, V, Comm>;
                    using Smoother = Dune::BlockPreconditioner<V, V, Comm, Dune::SeqSSOR<P,V,V>>;
                    auto lopP = std::make_shared<PressOp>(pmatrix_, *presComm_);
                    preconditionerForP_ = std::make_shared<
                        Dune::Amg::AMG<PressOp, V, Smoother, Comm>>(lopP, params, *presComm_);
                }
            }
            else
#endif
            {
                auto lopP = std::make_shared<PressureLinearOperator>(pmatrix_);
                preconditionerForP_ = std::make_shared<
                    Dune::Amg::AMG<PressureLinearOperator, V, Dumux::ParMTSSOR<P,V,V>>
                >(lopP, params);
            }
        }
        else
        {
            const std::size_t numIterations = pmatrix_.nonzeroes() == pmatrix_.N() ? 1 : 10;
#if HAVE_MPI
            if (presComm_)
            {
                const bool nonovP = (presComm_->category()
                                      == Dune::SolverCategory::nonoverlapping);
                auto seqJacP = std::make_shared<Dune::SeqJac<P,V,V>>(pmatrix_, numIterations, 1.0);
                if (nonovP)
                {
                    using PressJacobi = Dune::NonoverlappingBlockPreconditioner<Comm, Dune::SeqJac<P,V,V>>;
                    preconditionerForP_ = std::make_shared<PressJacobi>(seqJacP, *presComm_);
                }
                else
                {
                    using PressJacobi = Dune::BlockPreconditioner<V, V, Comm, Dune::SeqJac<P,V,V>>;
                    preconditionerForP_ = std::make_shared<PressJacobi>(seqJacP, *presComm_);
                }
            }
            else
#endif
            {
                using PressJacobi = Dumux::ParMTJac<P, V, V>;
                preconditionerForP_ = std::make_shared<PressJacobi>(pmatrix_, numIterations, 1.0);
            }
        }
    }

    template<class Sol, class Rhs>
    void applyPreconditionerForA_(Sol& sol, Rhs& rhs) const
    {
        preconditionerForA_->apply(sol, rhs);
    }

    template<class Sol, class Rhs>
    void applyPreconditionerForP_(Sol& sol, Rhs& rhs) const
    {
        preconditionerForP_->apply(sol, rhs);
    }

public:
    /*!
     * \brief Add a diagonal to the inverse of the pressure block
     * \note For a transient problem the pressure operator approximates the Schur complement of the storage part
     *       of the velocity block. The viscous part contributes a mass matrix scaled by the inverse viscosity,
     *       whose inverse is this diagonal; both inverses are summed (Cahouet & Chabard,
     *       Int. J. Numer. Meth. Fluids 8 (1988) 869-895).
     */
    void setPressureDiagonal(std::shared_ptr<const V> diagonal)
    { pressureDiagonal_ = std::move(diagonal); }

private:
    std::shared_ptr<const V> pressureDiagonal_;

    //! \brief The matrix we operate on.
    const M& matrix_;
    //! \brief The matrix we operate on.
    const P& pmatrix_;
    //! \brief The verbosity level
    const int verbosity_;

    std::shared_ptr<Dune::Preconditioner<U, U>> preconditionerForA_;
    std::shared_ptr<Dune::Preconditioner<V, V>> preconditionerForP_;
    std::shared_ptr<Dune::InverseOperator<U, U>> directSolverForA_;
    std::shared_ptr<Dune::InverseOperator<V, V>> directSolverForP_;

    const std::string paramGroup_;
    Mode mode_;

#if HAVE_MPI
    std::shared_ptr<Comm> velComm_;
    std::shared_ptr<Comm> presComm_;
#endif
};

#if HAVE_MPI

/*!
 * \brief Parallel linear operator for the Stokes saddle-point system
 *
 * The result is unique. On a non-overlapping decomposition ghost values propagate through
 * identity rows as in the NonoverlappingSchwarzOperator.
 */
template<class M, class X, class Y>
class ParallelStokesLinearOperator : public Dune::LinearOperator<X, Y>
{
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
public:
    using matrix_type = M;
    using domain_type = X;
    using range_type  = Y;
    using field_type  = typename X::field_type;

    ParallelStokesLinearOperator(const M& A,
                                  std::shared_ptr<Comm> vComm,
                                  std::shared_ptr<Comm> pComm,
                                  bool nonOverlapping)
    : A_(A), localOp_(A), vComm_(std::move(vComm)), pComm_(std::move(pComm))
    , nonOverlapping_(nonOverlapping)
    {
#if HAVE_MPI
        using namespace Dune::Indices;
        if (nonOverlapping_)
            velOp_ = std::make_shared<VelOp>(A_[_0][_0], *vComm_);
#endif
    }

    void apply(const X& x, Y& y) const override
    {
        using namespace Dune::Indices;
        if (nonOverlapping_)
        {
            applyNonoverlapping_(x, y);
        }
        else
        {
            // on an overlapping decomposition only the rows of owned dofs are complete
            localOp_.apply(x, y);
            makeUnique(*vComm_, y[_0]);
            makeUnique(*pComm_, y[_1]);
        }
    }

    void applyscaleadd(field_type alpha, const X& x, Y& y) const override
    {
        using namespace Dune::Indices;
        Y Ax(y); Ax = 0;
        apply(x, Ax);
        y[_0].axpy(alpha, Ax[_0]);
        y[_1].axpy(alpha, Ax[_1]);
    }

private:
    using A00t = std::decay_t<decltype(std::declval<M>()[Dune::Indices::_0][Dune::Indices::_0])>;
    using Ut   = std::decay_t<decltype(std::declval<X>()[Dune::Indices::_0])>;
    using VelOp = Dune::NonoverlappingSchwarzOperator<A00t, Ut, Ut, Comm>;

    //! The velocity block is summed over the processes and applied as NonoverlappingSchwarzOperator.
    //! The off-diagonal blocks stay additive: the rows of border dofs can reach dofs that only a
    //! neighbouring process assembles, so their partial products are summed instead.
    void applyNonoverlapping_(const X& x, Y& y) const
    {
        using namespace Dune::Indices;

        X xc(x);
        vComm_->copyOwnerToAll(xc[_0], xc[_0]);
        pComm_->copyOwnerToAll(xc[_1], xc[_1]);

        velOp_->apply(xc[_0], y[_0]);

        Ut bv(y[_0]); bv = 0.0;
        A_[_0][_1].umv(xc[_1], bv);
        zeroOverlap(*vComm_, bv);
        vComm_->addOwnerCopyToOwnerCopy(bv, bv);
        y[_0] += bv;

        // the pressure block is nonzero only in constrained rows, e.g. a pressure pin
        y[_1] = 0.0;
        A_[_1][_0].umv(xc[_0], y[_1]);
        A_[_1][_1].umv(xc[_1], y[_1]);
        zeroOverlap(*pComm_, y[_1]);
        pComm_->addOwnerCopyToOwnerCopy(y[_1], y[_1]);

        // the defect b - Ax the preconditioner receives has to be unique (Blatt & Bastian §4.3)
        makeUnique(*vComm_, y[_0]);
        makeUnique(*pComm_, y[_1]);
    }

public:

    Dune::SolverCategory::Category category() const override
    {
        return nonOverlapping_ ? Dune::SolverCategory::nonoverlapping
                               : Dune::SolverCategory::overlapping;
    }

private:
    const M& A_;
    Dumux::ParallelMultiTypeMatrixAdapter<M, X, Y> localOp_;
    std::shared_ptr<Comm> vComm_, pComm_;
    bool nonOverlapping_;
    std::shared_ptr<VelOp> velOp_;
};

/*!
 * \brief Parallel wrapper of the Stokes saddle-point preconditioner
 *
 * Applies the StokesPreconditioner and makes its result consistent.
 */
template<class SeqPrec, class X, class Y>
class ParallelStokesPreconditioner : public Dune::Preconditioner<X, Y>
{
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
public:
    using domain_type = X;
    using range_type  = Y;
    using field_type  = typename X::field_type;

    ParallelStokesPreconditioner(std::shared_ptr<SeqPrec> seqPrec,
                                  std::shared_ptr<Comm> vComm,
                                  std::shared_ptr<Comm> pComm,
                                  bool nonOverlapping)
    : seqPrec_(std::move(seqPrec))
    , vComm_(std::move(vComm)), pComm_(std::move(pComm))
    , nonOverlapping_(nonOverlapping)
    {}

    void pre(X& v, Y& d) override
    {
        seqPrec_->pre(v, d);
        communicateUpdate_(v);
    }

    void apply(X& v, const Y& d) override
    {
        seqPrec_->apply(v, d);
        communicateUpdate_(v);
    }

    void post(X& v) override
    {
        seqPrec_->post(v);
        communicateUpdate_(v);
    }

    Dune::SolverCategory::Category category() const override
    {
        return nonOverlapping_ ? Dune::SolverCategory::nonoverlapping
                               : Dune::SolverCategory::overlapping;
    }

private:
    void communicateUpdate_(X& v)
    {
        using namespace Dune::Indices;
        // the Krylov solver works with consistent vectors; dofs that are not shared,
        // such as bubble dofs, are not in the index set and need no communication
        vComm_->copyOwnerToAll(v[_0], v[_0]);
        pComm_->copyOwnerToAll(v[_1], v[_1]);
    }

    std::shared_ptr<SeqPrec> seqPrec_;
    std::shared_ptr<Comm> vComm_, pComm_;
    bool nonOverlapping_;
};

/*!
 * \brief Parallel scalar product for the velocity-pressure vector of the Stokes system
 *
 * Every dof contributes once, on the process that owns it.
 */
template<class X>
class ParallelStokesScalarProduct : public Dune::ScalarProduct<X>
{
    using Comm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
public:
    using domain_type  = X;
    using field_type   = typename X::field_type;
    using real_type    = typename Dune::FieldTraits<field_type>::real_type;

    ParallelStokesScalarProduct(std::shared_ptr<Comm> vComm,
                                 std::shared_ptr<Comm> pComm,
                                 bool nonOverlapping)
    : vComm_(std::move(vComm)), pComm_(std::move(pComm))
    , nonOverlapping_(nonOverlapping)
    {}

    field_type dot(const X& x, const X& y) const override
    {
        using namespace Dune::Indices;
        using Attr = Dune::OwnerOverlapCopyAttributeSet;

        // dofs that are not shared, such as bubble dofs, are not in the index set,
        // so sum all local entries and subtract those of shared dofs owned elsewhere
        field_type local = 0;

        for (std::size_t i = 0; i < x[_0].size(); ++i)
            local += Dune::dot(x[_0][i], y[_0][i]);
        for (const auto& pair : vComm_->indexSet())
            if (pair.local().attribute() != Attr::owner)
                local -= Dune::dot(x[_0][pair.local()], y[_0][pair.local()]);

        for (std::size_t i = 0; i < x[_1].size(); ++i)
            local += Dune::dot(x[_1][i], y[_1][i]);
        for (const auto& pair : pComm_->indexSet())
            if (pair.local().attribute() != Attr::owner)
                local -= Dune::dot(x[_1][pair.local()], y[_1][pair.local()]);

        return vComm_->communicator().sum(local);
    }

    real_type norm(const X& x) const override
    {
        using std::sqrt;
        return sqrt(std::real(dot(x, x)));
    }

    Dune::SolverCategory::Category category() const override
    {
        return nonOverlapping_ ? Dune::SolverCategory::nonoverlapping
                               : Dune::SolverCategory::overlapping;
    }


private:
    std::shared_ptr<Comm> vComm_, pComm_;
    bool nonOverlapping_;
};

#endif // HAVE_MPI

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Preconditioned iterative solver for the incompressible Stokes problem
 * \note Uses StokesPreconditioner as preconditioner (tailored to the incompressible Stokes problem)
 * \note Sequential and MPI-parallel (non-overlapping or overlapping decomposition)
 */
template<class Matrix, class Vector, class VelocityGG, class PressureGG>
class StokesSolver
: public LinearSolver
{
    using Preconditioner = Detail::StokesPreconditioner<Matrix, Vector, Vector>;
    using PressureVector = std::decay_t<decltype(std::declval<Vector>()[Dune::Indices::_1])>;

#if HAVE_MPI
    using VTraits = LinearSolverTraits<VelocityGG>;
    using PTraits = LinearSolverTraits<PressureGG>;
    using Comm    = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
    static constexpr bool canCommunicate_ = VTraits::canCommunicate && PTraits::canCommunicate;
#endif

public:
    //! The matrix type of the pressure operator of the preconditioner
    using PressureMatrix = typename Preconditioner::PressureLinearOperator::matrix_type;

    /*!
     * \brief Constructor
     * \param vGridGeometry grid geometry of the velocity discretization
     * \param pGridGeometry grid geometry of the pressure discretization
     * \param dirichletDofs a vector (same size and shape as right hand side) where the dirichlet dofs are marked with 1.0
     *                      and all other entries are 0.0.
     * \param paramGroup group prefix when looking up keys in the parameter tree
     */
    StokesSolver(std::shared_ptr<const VelocityGG> vGridGeometry,
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

        if (hasParamInGroup(this->paramGroup(), "Component.LiquidDynamicViscosity"))
            viscosity_ = getParamFromGroup<double>(this->paramGroup(), "Component.LiquidDynamicViscosity");
        else if (hasParamInGroup(this->paramGroup(), "Component.LiquidKinematicViscosity"))
            viscosity_ = getParamFromGroup<double>(this->paramGroup(), "Component.LiquidKinematicViscosity") * density_;
        else
        {
            const std::string group = this->paramGroup() == "" ? "Component" : this->paramGroup() + ".Component";
            DUNE_THROW(ParameterException, "Stokes solver requires the parameters"
                           << " LiquidDynamicViscosity or LiquidKinematicViscosity"
                           << " in parameter group [" << group << "]");
        }

        weight_     = getParamFromGroup<double>(this->paramGroup(), "LinearSolver.Preconditioner.MassMatrixWeight", 1.0);
        solverType_ = getParamFromGroup<std::string>(this->paramGroup(), "LinearSolver.Type", "gmres");

        scalarProduct_ = std::make_shared<Dune::ScalarProduct<Vector>>();
#if HAVE_MPI
        if constexpr (canCommunicate_)
        {
            if (vGridGeometry_->gridView().comm().size() > 1)
                initParallelInfrastructure_();
        }
#endif
    }

    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        auto bTmp = b;
        auto ATmp = A;

        return applyIterativeSolver_(ATmp, x, bTmp);
    }

    /*!
     * \brief Keep a copy of the matrix, and the operator and preconditioner built on it, for solve(x, b)
     * \note The preconditioner uses the pressure operator and diagonal set at the time of this call.
     */
    void setMatrix(const Matrix& A)
    {
        using namespace Dune::Indices;
        // the operator and preconditioner refer to the kept matrix
        keptComponents_ = {};
        keptMatrix_ = std::make_shared<Matrix>(A);
        eliminatedEntries_.reset();

        if (symmetrizeDirichlet_())
        {
            // the symmetrization moves the known Dirichlet values to the right-hand side; the entries
            // it eliminates are kept to do the same for every right-hand side passed to solve(x, b)
            auto unusedRhs = dirichletDofs_;
            symmetrizeConstraints(*keptMatrix_, unusedRhs, dirichletDofs_);
            eliminatedEntries_ = std::make_shared<Matrix>(A);
            *eliminatedEntries_ -= *keptMatrix_;
        }

        (*keptMatrix_)[_1] *= -1.0/density_;

#if HAVE_MPI
        if constexpr (canCommunicate_)
        {
            if (vGridGeometry_->gridView().comm().size() > 1 && isNonOverlapping_)
                prepareParallelMatrix_(*keptMatrix_);
        }
#endif

        keptComponents_ = makeSolverComponents_(*keptMatrix_);
    }

    //! \copydoc setMatrix(const Matrix&)
    void setMatrix(std::shared_ptr<Matrix> A)
    { setMatrix(*A); }

    //! Solve with the matrix kept by setMatrix
    bool solve(Vector& x, const Vector& b)
    {
        if (!keptMatrix_)
            DUNE_THROW(Dune::InvalidStateException, "Called solve(x, b) but no matrix has been set");

        using namespace Dune::Indices;
        auto bTmp = b;
        if (eliminatedEntries_)
        {
            auto dirichletValues = b;
            Dune::Hybrid::forEach(std::make_index_sequence<Vector::size()>{}, [&](auto i)
            {
                for (std::size_t k = 0; k < dirichletValues[i].size(); ++k)
                    for (std::size_t c = 0; c < dirichletValues[i][k].size(); ++c)
                        dirichletValues[i][k][c] *= dirichletDofs_[i][k][c];
            });
            eliminatedEntries_->mmv(dirichletValues, bTmp);
        }

        bTmp[_1] *= -1.0/density_;

#if HAVE_MPI
        if constexpr (canCommunicate_)
        {
            if (vGridGeometry_->gridView().comm().size() > 1)
            {
                prepareParallelRhs_(bTmp);
                const bool converged = runSolver_(keptComponents_.op, scalarProduct_, keptComponents_.preconditioner, x, bTmp);
                makeConsistent(x);
                return converged;
            }
        }
#endif
        return runSolver_(keptComponents_.op, scalarProduct_, keptComponents_.preconditioner, x, bTmp);
    }

    /*!
     * \brief Use the given matrix in place of the pressure mass matrix in the preconditioner
     * \note The matrix has to approximate the pressure Schur complement of the system. The viscosity-weighted
     *       mass matrix does so for a Stokes problem dominated by viscosity; a transient problem with a small
     *       time step is dominated by the storage term of the velocity block and calls for a Poisson-type
     *       operator weighted by the time step size and the inverse density.
     * \note On a non-overlapping decomposition the rows of border dofs have to be complete.
     */
    void setPressureMatrix(std::shared_ptr<const PressureMatrix> matrix)
    { pressureMatrix_ = std::move(matrix); }

    //! Set the diagonal added to the inverse of the pressure block (see StokesPreconditioner::setPressureDiagonal)
    void setPressureDiagonal(std::shared_ptr<const PressureVector> diagonal)
    { pressureDiagonal_ = std::move(diagonal); }

    //! Give the copies of shared degrees of freedom the values of their owners
    void makeConsistent(Vector& v) const
    {
#if HAVE_MPI
        using namespace Dune::Indices;
        if (vComm_ && pComm_)
        {
            vComm_->copyOwnerToAll(v[_0], v[_0]);
            pComm_->copyOwnerToAll(v[_1], v[_1]);
        }
#endif
    }

    //! Set the entries of shared degrees of freedom this process does not own to zero
    void makeUnique(Vector& v) const
    {
#if HAVE_MPI
        using namespace Dune::Indices;
        if (vComm_ && pComm_)
        {
            Detail::makeUnique(*vComm_, v[_0]);
            Detail::makeUnique(*pComm_, v[_1]);
        }
#endif
    }

    Scalar norm(const Vector& b) const
    {
#if HAVE_MPI
        // on a non-overlapping decomposition an assembled residual is additive, so the owned entries
        // are only complete after summing over the processes; on an overlapping one they are complete
        using namespace Dune::Indices;
        if (vComm_ && pComm_ && vGridGeometry_->gridView().comm().size() > 1)
        {
            Vector bc(b);
            if (isNonOverlapping_)
            {
                vComm_->addOwnerCopyToOwnerCopy(bc[_0], bc[_0]);
                pComm_->addOwnerCopyToOwnerCopy(bc[_1], bc[_1]);
            }
            else
            {
                Detail::makeUnique(*vComm_, bc[_0]);
                Detail::makeUnique(*pComm_, bc[_1]);
            }
            return scalarProduct_->norm(bc);
        }
#endif
        return scalarProduct_->norm(b);
    }

    std::string name() const
    {
        return "Block-preconditioned Stokes solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
#if HAVE_MPI
    void initParallelInfrastructure_()
    {
        // if either sub-problem is overlapping, use overlapping for both
        isNonOverlapping_ = VTraits::isNonOverlapping(vGridGeometry_->gridView())
                         && PTraits::isNonOverlapping(pGridGeometry_->gridView());

        const auto cat = isNonOverlapping_ ? Dune::SolverCategory::nonoverlapping
                                           : Dune::SolverCategory::overlapping;

        // the helpers keep a reference to the dof mapper, which the traits may return by value
        vDofMapper_ = std::make_shared<const typename VTraits::DofMapper>(VTraits::dofMapper(*vGridGeometry_));
        pDofMapper_ = std::make_shared<const typename PTraits::DofMapper>(PTraits::dofMapper(*pGridGeometry_));
        vParallelHelper_ = std::make_shared<ParallelISTLHelper<VTraits>>(vGridGeometry_->gridView(), *vDofMapper_);
        pParallelHelper_ = std::make_shared<ParallelISTLHelper<PTraits>>(pGridGeometry_->gridView(), *pDofMapper_);

        vComm_ = std::make_shared<Comm>(vGridGeometry_->gridView().comm(), cat);
        vParallelHelper_->createParallelIndexSet(*vComm_);

        pComm_ = std::make_shared<Comm>(pGridGeometry_->gridView().comm(), cat);
        pParallelHelper_->createParallelIndexSet(*pComm_);

        scalarProduct_ = std::make_shared<Detail::ParallelStokesScalarProduct<Vector>>(
            vComm_, pComm_, isNonOverlapping_);
    }
#endif

    bool symmetrizeDirichlet_() const
    { return getParamFromGroup<bool>(this->paramGroup(), "LinearSolver.SymmetrizeDirichlet", true); }

    bool applyIterativeSolver_(Matrix& A, Vector& x, Vector& b)
    {
        // make Dirichlet boundary conditions symmetric
        if (symmetrizeDirichlet_())
            symmetrizeConstraints(A, b, dirichletDofs_);

        // make Matrix symmetric on the block-scale
        using namespace Dune::Indices;
        A[_1] *= -1.0/density_;
        b[_1] *= -1.0/density_;

#if HAVE_MPI
        if constexpr (canCommunicate_)
        {
            if (vGridGeometry_->gridView().comm().size() > 1)
            {
                if (isNonOverlapping_)
                    prepareParallelMatrix_(A);
                prepareParallelRhs_(b);
                const auto components = makeSolverComponents_(A);
                const bool converged = runSolver_(components.op, scalarProduct_, components.preconditioner, x, b);
                makeConsistent(x);
                return converged;
            }
        }
#endif
        const auto components = makeSolverComponents_(A);
        return runSolver_(components.op, scalarProduct_, components.preconditioner, x, b);
    }

    //! What a solve needs besides the matrix
    struct SolverComponents
    {
        std::shared_ptr<Dune::LinearOperator<Vector, Vector>> op;
        std::shared_ptr<Dune::Preconditioner<Vector, Vector>> preconditioner;
        // the preconditioner refers to the matrix of the pressure operator
        std::shared_ptr<typename Preconditioner::PressureLinearOperator> pressureOperator;
    };

    SolverComponents makeSolverComponents_(const Matrix& A)
    {
        SolverComponents components;
        components.pressureOperator = makePressureLinearOperator_<typename Preconditioner::PressureLinearOperator>();
        auto assembledOperator = std::make_shared<Dumux::ParallelMultiTypeMatrixAdapter<Matrix, Vector, Vector>>(A);

#if HAVE_MPI
        if constexpr (canCommunicate_)
        {
            if (vGridGeometry_->gridView().comm().size() > 1)
            {
                // the preconditioner needs the assembled operator, the Krylov solver the parallel one
                components.op = std::make_shared<Detail::ParallelStokesLinearOperator<Matrix, Vector, Vector>>(
                    A, vComm_, pComm_, isNonOverlapping_);
                auto preconditioner = std::make_shared<Preconditioner>(
                    assembledOperator, components.pressureOperator, params_.sub("preconditioner"), vComm_, pComm_);
                preconditioner->setPressureDiagonal(pressureDiagonal_);
                components.preconditioner = std::make_shared<Detail::ParallelStokesPreconditioner<Preconditioner, Vector, Vector>>(
                    preconditioner, vComm_, pComm_, isNonOverlapping_);
                return components;
            }
        }
#endif

        auto preconditioner = std::make_shared<Preconditioner>(assembledOperator, components.pressureOperator, params_.sub("preconditioner"));
        preconditioner->setPressureDiagonal(pressureDiagonal_);
        components.op = assembledOperator;
        components.preconditioner = preconditioner;
        return components;
    }

#if HAVE_MPI
    void prepareParallelMatrix_(Matrix& A)
    {
        using namespace Dune::Indices;
        using GridView = typename VelocityGG::GridView;
        static constexpr int dim = GridView::dimension;
        static constexpr std::size_t numCodims = dim + 1;

        const auto& gv        = vGridGeometry_->gridView();
        const auto velCodims  = activeCodimsBitset_<VTraits, numCodims>();

        // the velocity block is summed over the processes: the smoother of the velocity AMG needs
        // complete diagonal entries; the off-diagonal blocks stay additive (see ParallelStokesLinearOperator)
        using VelBlock = std::decay_t<decltype(A[_0][_0])>;
        MultiCodimParallelMatrixHelper<VelBlock, GridView, typename VTraits::DofMapper, numCodims>
            velHelper(gv, *vDofMapper_, velCodims);
        velHelper.extendMatrix(A[_0][_0], [](auto){ return false; });
        velHelper.sumEntries(A[_0][_0]);
    }

    //! Bring the right-hand side into the unique representation the operator's result has
    void prepareParallelRhs_(Vector& b)
    {
        using namespace Dune::Indices;
        if (isNonOverlapping_)
        {
            using VelVec  = std::decay_t<decltype(b[_0])>;
            using PresVec = std::decay_t<decltype(b[_1])>;
            using VNonoverlapping = typename VTraits::template ParallelNonoverlapping<
                std::decay_t<decltype(std::declval<Matrix>()[_0][_0])>, VelVec>;
            using PNonoverlapping = typename PTraits::template ParallelNonoverlapping<
                std::decay_t<decltype(std::declval<Matrix>()[_1][_1])>, PresVec>;
            prepareVectorParallel<VTraits, VNonoverlapping>(b[_0], *vParallelHelper_);
            prepareVectorParallel<PTraits, PNonoverlapping>(b[_1], *pParallelHelper_);
        }

        Detail::makeUnique(*vComm_, b[_0]);
        Detail::makeUnique(*pComm_, b[_1]);
    }

    // Returns active codims as bitset<numCodims> from a LinearSolverTraits type.
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
        std::unique_ptr<Dune::InverseOperator<Vector, Vector>> solver;
        if (solverType_ == "minres")
            solver = std::make_unique<Dune::MINRESSolver<Vector>>(op, sp, prec, params_);
        else if (solverType_ == "bicgstab")
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

    template<class LinearOperator>
    std::shared_ptr<LinearOperator> makePressureLinearOperator_()
    {
        if (pressureMatrix_)
            return std::make_shared<LinearOperator>(pressureMatrix_);

        using M = typename LinearOperator::matrix_type;
        auto massMatrix = createMassMatrix_<M>();
        return std::make_shared<LinearOperator>(massMatrix);
    }

    template<class M, class DM = typename PressureGG::DiscretizationMethod>
    std::shared_ptr<M> createMassMatrix_()
    {
        auto massMatrix = std::make_shared<M>();
        massMatrix->setBuildMode(M::random);
        const auto numDofs = pGridGeometry_->numDofs();

        if constexpr (DM{} == DiscretizationMethods::cctpfa || DM{} == DiscretizationMethods::ccmpfa)
        {
            Dune::MatrixIndexSet pattern;
            pattern.resize(numDofs, numDofs);
            for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
                pattern.add(globalI, globalI);
            pattern.exportIdx(*massMatrix);

            const auto& gv = pGridGeometry_->gridView();
            auto fvGeometry = localView(*pGridGeometry_);
            for (const auto& element : elements(gv))
            {
                fvGeometry.bindElement(element);
                for (const auto& scv : scvs(fvGeometry))
                {
                    using Extrusion = Extrusion_t<PressureGG>;
                    const auto dofIndex = scv.dofIndex();
                    if (element.partitionType() == Dune::GhostEntity) // do not modify ghosts
                        (*massMatrix)[dofIndex][dofIndex][0][0] = 1.0;
                    else
                        (*massMatrix)[dofIndex][dofIndex][0][0] += weight_*Extrusion::volume(fvGeometry, scv)/(2.0*viscosity_);
                }
            }

            return massMatrix;
        }
        else if constexpr (DM{} == DiscretizationMethods::box || DM{} == DiscretizationMethods::fcdiamond)
        {
            getJacobianPattern<true>(*pGridGeometry_).exportIdx(*massMatrix);
            (*massMatrix) = 0.0;

            const auto& gv = pGridGeometry_->gridView();
            auto fvGeometry = localView(*pGridGeometry_);
            std::vector<Dune::FieldVector<double, 1>> values;

            for (const auto& element : elements(gv))
            {
                const auto geometry = element.geometry();
                const auto& localBasis = pGridGeometry_->feCache().get(element.type()).localBasis();
                const auto numLocalDofs = localBasis.size();
                values.resize(numLocalDofs);

                fvGeometry.bindElement(element);
                for (const auto& scvJ : scvs(fvGeometry))
                {
                    // Use mid-point rule (only works for linear ansatz functions)
                    const auto globalJ = scvJ.dofIndex();
                    const auto qWeightJ = PressureGG::Extrusion::volume(fvGeometry, scvJ);
                    const auto quadPos = geometry.local(scvJ.center());
                    localBasis.evaluateFunction(quadPos, values);

                    for (const auto& scvI : scvs(fvGeometry))
                    {
                        const auto valueIJ = values[scvI.localDofIndex()]*qWeightJ/(2.0*viscosity_);
                        (*massMatrix)[scvI.dofIndex()][globalJ][0][0] += valueIJ;
                    }
                }
            }

#if HAVE_MPI
            // on a non-overlapping decomposition the mass of a border dof is summed over the processes
            if (pGridGeometry_->gridView().comm().size() > 1 && isNonOverlapping_ && pParallelHelper_)
            {
                using PNonoverlapping = typename PTraits::template ParallelNonoverlapping<
                    M, std::decay_t<decltype((*massMatrix)[0])>>;
                prepareMatrixParallel<PTraits, PNonoverlapping>(*massMatrix, *pParallelHelper_);
            }
#endif

            return massMatrix;
        }

        DUNE_THROW(Dune::NotImplemented, "Mass matrix for discretization method not implemented");
    }

    double density_, viscosity_, weight_;
    Dune::InverseOperatorResult result_;
    Dune::ParameterTree params_;
    std::shared_ptr<const VelocityGG> vGridGeometry_;
    std::shared_ptr<const PressureGG> pGridGeometry_;
    const Vector& dirichletDofs_;
    std::string solverType_;

    std::shared_ptr<Dune::ScalarProduct<Vector>> scalarProduct_;

    std::shared_ptr<const PressureMatrix> pressureMatrix_;
    std::shared_ptr<const PressureVector> pressureDiagonal_;
    std::shared_ptr<Matrix> keptMatrix_;
    std::shared_ptr<Matrix> eliminatedEntries_;
    SolverComponents keptComponents_;

#if HAVE_MPI
    std::shared_ptr<const typename VTraits::DofMapper> vDofMapper_;
    std::shared_ptr<const typename PTraits::DofMapper> pDofMapper_;
    std::shared_ptr<ParallelISTLHelper<VTraits>> vParallelHelper_;
    std::shared_ptr<ParallelISTLHelper<PTraits>> pParallelHelper_;
    std::shared_ptr<Comm> vComm_, pComm_;
    bool isNonOverlapping_ = true;
#endif
};

} // end namespace Dumux

#endif
