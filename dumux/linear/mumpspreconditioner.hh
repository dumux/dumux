// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief A sequential preconditioner doing an exact local solve with MUMPS
 *
 * \ref Dumux::SeqMumps is a sequential preconditioner that factorizes the
 * (rank-)local matrix block once with MUMPS (on MPI_COMM_SELF) and applies the
 * exact inverse of that block in every preconditioner step. On its own this is a
 * direct solve. Its purpose is to be used as the local subdomain solve of a
 * parallel domain-decomposition preconditioner: when dune-istl wraps a sequential
 * preconditioner for an overlapping parallel operator it uses
 * Dune::BlockPreconditioner, whose apply() does the local solve followed by
 * copyOwnerToAll(). Combined with the exact local solve provided here this yields
 * a Restricted Additive Schwarz (RAS) preconditioner with exact subdomain solves.
 *
 * Select it through the ISTL solver factory (Dumux::IstlSolverFactoryBackend) by
 * setting LinearSolver.Preconditioner.Type = seqmumps. In an overlapping parallel
 * run (e.g. cell-centered schemes with Grid.Overlap >= 1) the factory automatically
 * wraps it into the RAS block preconditioner; sequentially it is an exact solve.
 *
 * Unlike Dumux::DirectSolverMumps (a full distributed parallel direct solver) this
 * class only ever operates on the local matrix using a private MPI_COMM_SELF
 * communicator, so several ranks each factorize their own subdomain independently.
 */
#ifndef DUMUX_LINEAR_MUMPS_PRECONDITIONER_HH
#define DUMUX_LINEAR_MUMPS_PRECONDITIONER_HH

#include <vector>
#include <limits>
#include <memory>
#include <cstring>

#include <dune/common/exceptions.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/simd/simd.hh>

#include <dune/istl/preconditioner.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvercategory.hh>

#include <dumux/linear/istlsolverregistry.hh>

#if DUMUX_HAVE_MUMPS

#include <mpi.h>
#include <dmumps_c.h>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Sequential preconditioner performing an exact local solve with MUMPS
 *
 * \tparam M The (local) matrix type to operate on (a Dune::BCRSMatrix).
 * \tparam X Type of the update.
 * \tparam Y Type of the defect.
 * \tparam l Preconditioner block level (unused, exact solve; kept for the factory interface).
 */
template<class M, class X, class Y, int l = 1>
class SeqMumps : public Dune::Preconditioner<X, Y>
{
    using BlockType = typename M::block_type;
    static constexpr int blockRows = BlockType::rows;
    static constexpr int blockCols = BlockType::cols;

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
     * \brief Constructor from an assembled linear operator (used by the solver factory).
     */
    SeqMumps(const std::shared_ptr<const Dune::AssembledLinearOperator<M, X, Y>>& op,
             const Dune::ParameterTree& config)
    : SeqMumps(op->getmat(), config.get<int>("verbosity", 0))
    {}

    /*!
     * \brief Constructor from the local matrix. Factorizes the matrix once.
     */
    SeqMumps(const M& A, int verbosity = 0)
    : verbosity_(verbosity)
    {
        initMumps_();
        factorize_(A);
    }

    ~SeqMumps()
    {
        if (initialized_)
        {
            id_.job = -2; // terminate the MUMPS instance
            dmumps_c(&id_);
        }
    }

    // Non-copyable: owns a MUMPS instance and raw pointers into our buffers.
    SeqMumps(const SeqMumps&) = delete;
    SeqMumps& operator=(const SeqMumps&) = delete;

    /*!
     * \brief Re-factorize after the matrix values changed but the sparsity pattern did not.
     *
     * Reuses the cached symbolic analysis (MUMPS JOB=2), which is the expensive,
     * pattern-only step. Falls back to a full analysis + factorization if the structure
     * changed (different number of dofs or nonzeros). Use this to reuse one factorization
     * setup across e.g. Newton steps, where the Jacobian values change but the stencil does not.
     */
    void refactorize(const M& A)
    {
        if (!refillValues_(A))
        {
            factorize_(A); // structure changed -> rebuild triplets and re-analyze
            return;
        }
        id_.a = a_.data();
        id_.job = 2; // numeric factorization reusing the stored analysis
        dmumps_c(&id_);
        checkError_("numeric factorization");
    }

    void pre(X&, Y&) override {}

    /*!
     * \brief Apply the exact inverse of the local block: v = A^{-1} d.
     */
    void apply(X& v, const Y& d) override
    {
        // pack the defect into the dense centralized right-hand side (scalar order)
        std::size_t k = 0;
        for (std::size_t i = 0; i < d.size(); ++i)
            for (int ki = 0; ki < blockRows; ++ki)
                rhs_[k++] = d[i][ki];

        id_.nrhs = 1;
        id_.lrhs = id_.n;
        id_.rhs = rhs_.data();
        id_.job = 3; // solve (reusing the stored factorization)
        dmumps_c(&id_);
        checkError_("solve");

        // unpack the solution into the update
        k = 0;
        for (std::size_t i = 0; i < v.size(); ++i)
            for (int ki = 0; ki < blockRows; ++ki)
                v[i][ki] = rhs_[k++];
    }

    void post(X&) override {}

    //! Category of the preconditioner (see SolverCategory::Category)
    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:
    /*!
     * \brief Initialize a MUMPS instance on a private MPI_COMM_SELF communicator.
     *
     * sym=0: unsymmetric. par=1: the (single) host participates. ICNTL(18)=0:
     * centralized assembled input (we supply the whole local matrix on this rank).
     */
    void initMumps_()
    {
        std::memset(&id_, 0, sizeof(id_));
        id_.comm_fortran = static_cast<int>(MPI_Comm_c2f(MPI_COMM_SELF));
        id_.par = 1;
        id_.sym = 0;
        id_.job = -1;
        dmumps_c(&id_);
        initialized_ = true;
        checkError_("initialization");

        icntl_(1) = -1; // suppress error output stream
        icntl_(2) = -1; // suppress diagnostic/warning output stream
        icntl_(3) = -1; // suppress global information output stream
        icntl_(4) =  1; // verbosity: errors only
        icntl_(5) =  0; // assembled matrix (not elemental)
        icntl_(18) = 0; // centralized assembled input on the host
        icntl_(20) = 0; // dense, centralized RHS
        icntl_(21) = 0; // centralized solution
    }

    /*!
     * \brief Build the centralized triplet input and perform analysis + factorization.
     */
    void factorize_(const M& A)
    {
        const std::size_t n = A.N();
        nScalar_ = n * static_cast<std::size_t>(blockRows);

        std::size_t nnz = 0;
        for (auto row = A.begin(); row != A.end(); ++row)
            nnz += row->size() * static_cast<std::size_t>(blockRows) * blockCols;

        irn_.resize(nnz);
        jcn_.resize(nnz);
        a_.resize(nnz);
        rhs_.resize(nScalar_);

        std::size_t k = 0;
        for (auto row = A.begin(); row != A.end(); ++row)
        {
            const std::size_t i = row.index();
            for (int ki = 0; ki < blockRows; ++ki)
            {
                const std::size_t gRow = i * blockRows + ki;
                for (auto col = row->begin(); col != row->end(); ++col)
                {
                    const std::size_t j = col.index();
                    for (int kj = 0; kj < blockCols; ++kj)
                    {
                        irn_[k] = toMumpsInt_(gRow + 1);            // 1-based (Fortran)
                        jcn_[k] = toMumpsInt_(j * blockCols + kj + 1); // 1-based
                        a_[k] = (*col)[ki][kj];
                        ++k;
                    }
                }
            }
        }

        id_.n = toMumpsInt_(nScalar_);
        id_.nnz = static_cast<MUMPS_INT8>(a_.size());
        id_.nz = toMumpsInt_(a_.size()); // legacy field for older MUMPS
        id_.irn = irn_.data();
        id_.jcn = jcn_.data();
        id_.a = a_.data();

        id_.job = 4; // combined symbolic analysis (1) and numeric factorization (2)
        dmumps_c(&id_);
        checkError_("analysis and factorization");
    }

    /*!
     * \brief Refill the value array a_ from A, keeping the existing triplet structure.
     *
     * Returns false (without modifying a_) if the structure changed (different number of
     * scalar dofs or nonzeros), signalling that a full re-analysis is required. The iteration
     * order matches factorize_, so a_ stays consistent with irn_/jcn_.
     */
    bool refillValues_(const M& A)
    {
        if (A.N() * static_cast<std::size_t>(blockRows) != nScalar_)
            return false;

        std::size_t nnz = 0;
        for (auto row = A.begin(); row != A.end(); ++row)
            nnz += row->size() * static_cast<std::size_t>(blockRows) * blockCols;
        if (nnz != a_.size())
            return false;

        std::size_t k = 0;
        for (auto row = A.begin(); row != A.end(); ++row)
            for (int ki = 0; ki < blockRows; ++ki)
                for (auto col = row->begin(); col != row->end(); ++col)
                    for (int kj = 0; kj < blockCols; ++kj)
                        a_[k++] = (*col)[ki][kj];
        return true;
    }

    // Accessor for MUMPS control parameters using the manual's 1-based ICNTL(i) numbering.
    int& icntl_(int i) { return id_.icntl[i - 1]; }
    const int& infog_(int i) const { return id_.infog[i - 1]; }

    static MUMPS_INT toMumpsInt_(std::size_t v)
    {
        if (v > static_cast<std::size_t>(std::numeric_limits<MUMPS_INT>::max()))
            DUNE_THROW(Dune::RangeError, "Index " << v << " exceeds MUMPS_INT range; "
                       "rebuild MUMPS with 64-bit integers (-DINTSIZE64)");
        return static_cast<MUMPS_INT>(v);
    }

    void checkError_(const char* phase) const
    {
        if (infog_(1) < 0)
            DUNE_THROW(Dune::Exception, "SeqMumps " << phase << " failed: INFOG(1)="
                       << infog_(1) << ", INFOG(2)=" << infog_(2));
    }

    int verbosity_;
    DMUMPS_STRUC_C id_;
    bool initialized_ = false;
    std::size_t nScalar_ = 0;

    std::vector<MUMPS_INT> irn_;
    std::vector<MUMPS_INT> jcn_;
    std::vector<double> a_;
    std::vector<double> rhs_; // dense centralized RHS / solution
};

DUMUX_REGISTER_PRECONDITIONER("seqmumps", Dune::PreconditionerTag, Dune::defaultPreconditionerBlockLevelCreator<Dumux::SeqMumps>());

} // end namespace Dumux

#endif // DUMUX_HAVE_MUMPS

#endif
