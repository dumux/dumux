// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Parallel block-diagonal (Restricted Additive Schwarz) preconditioned BiCGSTAB
 *        with an exact local MUMPS solve per rank, for MultiType (multidomain) systems.
 *
 * This is the iterative, domain-decomposition counterpart of \ref Dumux::DirectSolverMumps
 * for coupled MultiTypeBlockMatrix systems. Instead of one global parallel factorization it
 * keeps the work local: the coupled 2x2 system is flattened to a scalar matrix (per rank, the
 * local part), each rank factorizes its own local block once per solve with sequential MUMPS
 * (\ref Dumux::SeqMumps on MPI_COMM_SELF), and a parallel BiCGSTAB iterates to the true global
 * solution. The preconditioner is Restricted Additive Schwarz (RAS): apply the local exact solve
 * and then keep only the owner's value of each shared dof (copyOwnerToAll).
 *
 * Why a hand-rolled Krylov: DuMux's generic ISTL parallel solvers fall back to sequential for
 * MultiTypeBlockMatrix, and a pure overlapping ISTL operator would be incorrect for subdomains
 * whose owner rows are incomplete after local assembly (e.g. a low-dimensional network
 * distributed by the bulk partition, which is effectively non-overlapping). We therefore reuse
 * the exact per-subdomain consistency machinery of DirectSolverMumps:
 *  - rhsPreprocess: sum partial border-dof contributions onto the owner (non-overlapping rows);
 *  - ghostBroadcast: propagate the owner's value to all copies (copyOwnerToAll).
 * The matrix-vector product and the (owned-only) inner products built from these closures define
 * a correct parallel Krylov method regardless of how good the local block preconditioner is.
 *
 * Requires MPI. Constructed with the tuple of subdomain grid geometries (same order as the
 * MultiTypeBlockVector subvectors), exactly like the multidomain DirectSolverMumps constructor.
 */
#ifndef DUMUX_LINEAR_MUMPS_RAS_SOLVER_HH
#define DUMUX_LINEAR_MUMPS_RAS_SOLVER_HH

#include <vector>
#include <tuple>
#include <memory>
#include <functional>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/preconditioners.hh> // Dune::SeqILU
#include <dune/istl/solver.hh> // Dune::InverseOperatorResult

#include <dumux/common/parameters.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/matrixconverter.hh>
#include <dumux/linear/mumpspreconditioner.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

#if DUMUX_HAVE_MUMPS

#include <mpi.h>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief RAS-preconditioned BiCGSTAB with exact local MUMPS block solves (multidomain).
 */
template<class LSTraits, class LATraits>
class MumpsRASSolver : public LinearSolver
{
    using Matrix = typename LATraits::Matrix;   // MultiTypeBlockMatrix
    using XVector = typename LATraits::Vector;  // MultiTypeBlockVector
    using BVector = typename LATraits::Vector;

    using Scalar = double;

    static constexpr bool isMultiType = isMultiTypeBlockVector<XVector>::value;
    static_assert(isMultiType, "MumpsRASSolver currently supports MultiType (multidomain) systems only");

    // Flattened scalar types produced by the converters (1 scalar dof per block).
    using ScalarBlockVector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;
    using ScalarMatrix = std::decay_t<decltype(
        MatrixConverter<Matrix>::multiTypeToBCRSMatrix(std::declval<Matrix>()))>;
    // The local subdomain solve is a Dune::Preconditioner on the flattened scalar block. It is
    // either an exact MUMPS factorization (SeqMumps) or an inexact ILU(k) (Dune::SeqILU).
    using MumpsLocal = SeqMumps<ScalarMatrix, ScalarBlockVector, ScalarBlockVector>;
    using LocalPreconditioner = Dune::Preconditioner<ScalarBlockVector, ScalarBlockVector>;

    // Per-subdomain parallel bookkeeping, in MultiTypeBlockVector subvector order. Mirrors the
    // closures of DirectSolverMumps and operates on the subdomain's contiguous scalar slice.
    struct SubdomainInfo
    {
        std::size_t scalarBS = 1;
        std::size_t nBlockLocal = 0;
        std::size_t scalarOffset = 0;
        bool rowsComplete = true;
        std::function<void(ScalarBlockVector&)> rhsPreprocess;  // sum partial border rows onto owner
        std::function<void(ScalarBlockVector&)> ghostBroadcast; // copyOwnerToAll (caller zeroes non-owned)
    };

public:
    struct SolverResult : public Dune::InverseOperatorResult
    {
        SolverResult() = default;
        SolverResult(const Dune::InverseOperatorResult& o) : InverseOperatorResult(o) {}
        operator bool() const { return this->converged; }
    };

    //! Sequential constructor (single rank). The local block is the whole system (exact solve).
    explicit MumpsRASSolver(const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
    , isParallel_(false)
    {
        mpiComm_ = MPI_COMM_SELF;
        MPI_Comm_rank(mpiComm_, &rank_);
        readParams_();
    }

    /*!
     * \brief Parallel multidomain constructor. \a ggTuple holds the subdomain grid geometries
     *        (or (shared) pointers to them), in the same order as the MultiTypeBlockVector.
     */
    template<class GGTuple>
    explicit MumpsRASSolver(const GGTuple& ggTuple, const std::string& paramGroup = "")
    requires (requires { std::tuple_size<GGTuple>::value; })
    : LinearSolver(paramGroup)
    , isParallel_(Detail::mumpsDeref(std::get<0>(ggTuple)).gridView().comm().size() > 1)
    {
        mpiComm_ = static_cast<MPI_Comm>(Detail::mumpsDeref(std::get<0>(ggTuple)).gridView().comm());
        MPI_Comm_rank(mpiComm_, &rank_);
        readParams_();
        setupLayout_(ggTuple);
    }

    SolverResult solve(const Matrix& A, XVector& x, const BVector& b)
    {
        // Flatten the coupled system to the (rank-local) scalar matrix and vectors.
        auto Ascalar = MatrixConverter<Matrix>::multiTypeToBCRSMatrix(A);
        auto bscalar = VectorConverter<BVector>::multiTypeToBlockVector(b);

        if (isOwned_.empty()) // sequential / lazy single-rank setup
            initSeqLayout_(Ascalar.N());

        // Prepare the local subdomain preconditioner (exact MUMPS or inexact ILU on this rank's
        // block). For MUMPS the symbolic analysis is cached across Newton steps (constant pattern).
        prepareLocal_(Ascalar);

        // Optional coarse-grid (two-level Schwarz) correction: assemble the small Galerkin coarse
        // operator A0 = R0 A R0^T for the current Jacobian.
        if (twoLevel_)
            buildCoarseOperator_(Ascalar);

        // Right-hand side: assemble partial border contributions onto the owner, then make it a
        // consistent global vector (every copy holds the owner's value).
        ScalarBlockVector bb = bscalar;
        sumPartialsToOwner_(bb);
        copyOwnerToAll_(bb);

        ScalarBlockVector xx(bb.size());
        xx = 0.0;

        auto result = bicgstab_(Ascalar, xx, bb);

        VectorConverter<XVector>::retrieveValues(x, xx);
        return result;
    }

    Scalar norm(const XVector& x) const
    {
        auto v = VectorConverter<XVector>::multiTypeToBlockVector(x);
        if (!isParallel_ || isOwned_.empty())
            return v.two_norm();

        // Parallel: assemble partial border contributions, then global norm over owned dofs.
        for (const auto& sub : subs_)
            if (sub.rhsPreprocess)
                sub.rhsPreprocess(v);

        double localSumSq = 0.0;
        for (std::size_t i = 0; i < v.size(); ++i)
            if (isOwned_[i])
                localSumSq += v[i][0] * v[i][0];

        double globalSumSq = 0.0;
        MPI_Allreduce(&localSumSq, &globalSumSq, 1, MPI_DOUBLE, MPI_SUM, mpiComm_);
        using std::sqrt;
        return sqrt(globalSumSq);
    }

    std::string name() const
    {
        return std::string("DuMux ") + (twoLevel_ ? "two-level " : "") + "RAS BiCGSTAB ("
             + (localSolverType_ == "ilu" ? "ILU" : "MUMPS") + " local solve)";
    }

private:
    // ---- parallel layout -------------------------------------------------------------------

    template<class GGTuple>
    void setupLayout_(const GGTuple& ggTuple)
    {
        constexpr std::size_t numSub = std::tuple_size_v<GGTuple>;
        static_assert(numSub == std::tuple_size_v<XVector>,
                      "Number of grid geometries must match the number of subdomains");
        subs_.resize(numSub);

        std::size_t scalarLocalOffset = 0;
        using Dune::Hybrid::forEach;
        forEach(std::make_index_sequence<numSub>{}, [&](auto I)
        {
            constexpr std::size_t i = decltype(I)::value;
            const auto& gg = Detail::mumpsDeref(std::get<i>(ggTuple));
            using SubGG = std::decay_t<decltype(gg)>;
            using SubTraits = LinearSolverTraits<SubGG>;
            using SubBlockVec = std::tuple_element_t<i, XVector>;
            constexpr std::size_t scalarBS = SubBlockVec::block_type::size();
            using GV = typename SubGG::GridView;
            constexpr bool canComm = subCanComm_<SubTraits, GV>();

            if (isParallel_ && !canComm)
                DUNE_THROW(Dune::NotImplemented,
                    "MumpsRASSolver: subdomain " << i << " has a DOF mapper that cannot drive the "
                    "generic grid communication; parallel RAS is not supported for it.");

            const auto& mapper = subDofMapper_<SubTraits>(gg);
            const std::size_t nBlockLocal = mapper.size();

            std::vector<bool> owned(nBlockLocal, true);
            bool rowsComplete = true;
            if (isParallel_)
            {
                if constexpr (canComm)
                {
                    ParallelISTLHelper<SubTraits> pHelper(gg.gridView(), mapper);
                    for (std::size_t b = 0; b < nBlockLocal; ++b)
                        owned[b] = pHelper.isOwned(b);

                    // Owner border rows complete after local (interior) assembly? Same rule as
                    // DirectSolverMumps: overlapping -> complete; non-overlapping -> incomplete for
                    // multi-codim / vertex-based schemes, complete otherwise (e.g. cell-centered).
                    const bool nonOverlapping = SubTraits::isNonOverlapping(gg.gridView());
                    constexpr bool hasMultiCodim = requires { SubTraits::dofCodims; };
                    constexpr bool vertexBased = (SubTraits::dofCodim == static_cast<int>(GV::dimension));
                    rowsComplete = nonOverlapping ? !(hasMultiCodim || vertexBased) : true;
                }
            }

            // Coarse-space class of each scalar dof: one constant basis function per subdomain
            // (field) and per equation component (used by the two-level correction).
            for (std::size_t b = 0; b < nBlockLocal; ++b)
                for (std::size_t k = 0; k < scalarBS; ++k)
                {
                    isOwned_.push_back(owned[b]);
                    coarseClass_.push_back(nLocalCoarse_ + k);
                }

            subs_[i].scalarBS = scalarBS;
            subs_[i].nBlockLocal = nBlockLocal;
            subs_[i].scalarOffset = scalarLocalOffset;
            subs_[i].rowsComplete = rowsComplete;
            if (isParallel_)
                if constexpr (canComm)
                    setupSubdomainClosures_<i>(gg);

            scalarLocalOffset += scalarBS * nBlockLocal;
            nLocalCoarse_ += scalarBS;
        });

        if (twoLevel_)
            setupCoarseSpace_();
    }

    void initSeqLayout_(std::size_t nScalar)
    {
        isOwned_.assign(nScalar, true);
        // single coarse class for the no-grid sequential path (two-level is a parallel feature)
        nLocalCoarse_ = 1;
        coarseClass_.assign(nScalar, 0);
        if (twoLevel_)
            setupCoarseSpace_();
    }

    // Build the rhsPreprocess / ghostBroadcast closures for subdomain i (copied from
    // DirectSolverMumps): they unpack the subdomain's scalar slice into a typed block vector,
    // run the Dune grid communication, and pack the result back.
    template<std::size_t i, class SubGG>
    void setupSubdomainClosures_(const SubGG& gg)
    {
        using SubTraits = LinearSolverTraits<SubGG>;
        using SubBlockVec = std::tuple_element_t<i, XVector>;
        using GV = typename SubGG::GridView;
        using Mapper = typename SubTraits::DofMapper;
        constexpr std::size_t scalarBS = SubBlockVec::block_type::size();

        const GV gv = gg.gridView();
        auto mapper = std::make_shared<Mapper>(subDofMapper_<SubTraits>(gg));
        const std::size_t offset = subs_[i].scalarOffset;
        const std::size_t nBlock = subs_[i].nBlockLocal;

        auto unpack = [offset, nBlock](const ScalarBlockVector& v, SubBlockVec& block) {
            block.resize(nBlock);
            for (std::size_t b = 0; b < nBlock; ++b)
                for (std::size_t k = 0; k < scalarBS; ++k)
                    block[b][k] = v[offset + scalarBS * b + k];
        };
        auto pack = [offset, nBlock](const SubBlockVec& block, ScalarBlockVector& v) {
            for (std::size_t b = 0; b < nBlock; ++b)
                for (std::size_t k = 0; k < scalarBS; ++k)
                    v[offset + scalarBS * b + k] = block[b][k];
        };

        // copyOwnerToAll: sum each dof's copies to the owner's value (callers zero non-owned first).
        subs_[i].ghostBroadcast = [gv, mapper, unpack, pack](ScalarBlockVector& v) {
            SubBlockVec block;
            unpack(v, block);
            if constexpr (requires { SubTraits::dofCodims; }) {
                MultiCodimVectorCommDataHandleSum<Mapper, SubBlockVec, GV::dimension>
                    handle(*mapper, block, SubTraits::dofCodims);
                gv.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
            } else {
                VectorCommDataHandleSum<Mapper, SubBlockVec, SubTraits::dofCodim>
                    handle(*mapper, block);
                gv.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
            }
            pack(block, v);
        };

        // rhsPreprocess: only for non-overlapping decompositions with incomplete owner rows.
        if (!subs_[i].rowsComplete && SubTraits::isNonOverlapping(gv)) {
            subs_[i].rhsPreprocess = [gv, mapper, unpack, pack](ScalarBlockVector& v) {
                SubBlockVec block;
                unpack(v, block);
                ParallelVectorHelper<GV, Mapper, SubTraits::dofCodim> vecHelper(gv, *mapper);
                if constexpr (requires { SubTraits::dofCodims; })
                    vecHelper.makeNonOverlappingConsistent(block, SubTraits::dofCodims);
                else
                    vecHelper.makeNonOverlappingConsistent(block);
                pack(block, v);
            };
        }
    }

    template<class SubTraits, class SubGG>
    static decltype(auto) subDofMapper_(const SubGG& gg)
    {
        if constexpr (requires { gg.dofMapper(); })
            return (gg.dofMapper());
        else
            return SubTraits::dofMapper(gg);
    }

    template<class SubTraits, class GV>
    static constexpr bool subCanComm_()
    {
        using Mapper = typename SubTraits::DofMapper;
        using Element = typename GV::template Codim<0>::Entity;
        return SubTraits::canCommunicate
            && requires (const Mapper& m, const Element& e) { m.indices(e); };
    }

    // ---- parallel vector operations on the flattened scalar vector -------------------------

    //! Sum partial border-row contributions onto the owner (non-overlapping subdomains).
    void sumPartialsToOwner_(ScalarBlockVector& v) const
    {
        if (!isParallel_) return;
        for (const auto& sub : subs_)
            if (sub.rhsPreprocess)
                sub.rhsPreprocess(v);
    }

    //! Make a vector consistent: each copy of a dof takes the owner's value.
    void copyOwnerToAll_(ScalarBlockVector& v) const
    {
        if (!isParallel_) return;
        for (std::size_t i = 0; i < v.size(); ++i)
            if (!isOwned_[i])
                v[i] = 0.0;
        for (const auto& sub : subs_)
            if (sub.ghostBroadcast)
                sub.ghostBroadcast(v);
    }

    //! Parallel matrix-vector product y = A x for a consistent x; returns a consistent y.
    void applyA_(const ScalarMatrix& A, const ScalarBlockVector& x, ScalarBlockVector& y) const
    {
        y.resize(x.size());
        y = 0.0;
        A.mv(x, y);             // local product (owner rows partial for non-overlapping subdomains)
        sumPartialsToOwner_(y); // assemble true owned values
        copyOwnerToAll_(y);     // make consistent
    }

    //! (Two-level restricted) additive Schwarz preconditioner.
    //! z = R^T (local solve) r  [+ R0^T A0^{-1} R0 r if the coarse correction is enabled].
    void applyM_(const ScalarBlockVector& r, ScalarBlockVector& z) const
    {
        z.resize(r.size());
        z = 0.0;
        localPrec_->apply(z, r); // local block solve (exact MUMPS or inexact ILU)
        copyOwnerToAll_(z);      // RAS restriction: keep the owner's value of each shared dof

        if (twoLevel_)
        {
            ScalarBlockVector zc;
            applyCoarse_(r, zc); // coarse-grid correction (additive)
            for (std::size_t i = 0; i < z.size(); ++i)
                z[i] += zc[i];
        }
    }

    //! Global inner product of two consistent vectors (counted once per dof, at its owner).
    Scalar dot_(const ScalarBlockVector& a, const ScalarBlockVector& b) const
    {
        double local = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i)
            if (isOwned_[i])
                local += a[i][0] * b[i][0];
        if (!isParallel_)
            return local;
        double global = 0.0;
        MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, mpiComm_);
        return global;
    }

    Scalar norm_(const ScalarBlockVector& a) const { using std::sqrt; return sqrt(dot_(a, a)); }

    // ---- preconditioned BiCGSTAB -----------------------------------------------------------

    SolverResult bicgstab_(const ScalarMatrix& A,
                           ScalarBlockVector& x, const ScalarBlockVector& b)
    {
        SolverResult result;
        const std::size_t maxIt = static_cast<std::size_t>(this->maxIter());
        const double reduction = this->residReduction();
        const bool verbose = this->verbosity() > 0 && rank_ == 0;

        ScalarBlockVector r(b), rHat(b.size()), p(b.size()), v(b.size()),
                          s(b.size()), t(b.size()), pHat(b.size()), sHat(b.size()), Ax(b.size());

        applyA_(A, x, Ax);          // x is 0 initially -> Ax = 0, but keep it general
        for (std::size_t i = 0; i < r.size(); ++i) r[i] -= Ax[i];

        rHat = r;
        p = 0.0; v = 0.0;
        double rho = 1.0, alpha = 1.0, omega = 1.0;

        const double bnorm = norm_(b);
        const double defScale = (bnorm > 0.0) ? bnorm : 1.0;
        double resNorm = norm_(r);

        if (verbose)
            std::cout << "=== MumpsRAS BiCGSTAB\n    Iter        Defect          Reduction\n"
                      << "       0   " << resNorm << "   " << resNorm/defScale << std::endl;

        result.converged = (resNorm <= reduction * defScale);
        std::size_t it = 0;
        for (; it < maxIt && !result.converged; ++it)
        {
            const double rhoNew = dot_(rHat, r);
            if (std::abs(rhoNew) < 1e-300) break; // breakdown
            const double beta = (rhoNew / rho) * (alpha / omega);
            for (std::size_t i = 0; i < p.size(); ++i)
                p[i] = r[i] + beta * (p[i] - omega * v[i]);

            applyM_(p, pHat);
            applyA_(A, pHat, v);
            const double rHatV = dot_(rHat, v);
            if (std::abs(rHatV) < 1e-300) break;
            alpha = rhoNew / rHatV;

            for (std::size_t i = 0; i < s.size(); ++i)
                s[i] = r[i] - alpha * v[i];

            const double sNorm = norm_(s);
            if (sNorm <= reduction * defScale)
            {
                for (std::size_t i = 0; i < x.size(); ++i)
                    x[i] += alpha * pHat[i];
                resNorm = sNorm;
                result.converged = true;
                ++it;
                break;
            }

            applyM_(s, sHat);
            applyA_(A, sHat, t);
            const double tt = dot_(t, t);
            if (std::abs(tt) < 1e-300) break;
            omega = dot_(t, s) / tt;

            for (std::size_t i = 0; i < x.size(); ++i)
                x[i] += alpha * pHat[i] + omega * sHat[i];
            for (std::size_t i = 0; i < r.size(); ++i)
                r[i] = s[i] - omega * t[i];

            rho = rhoNew;
            resNorm = norm_(r);
            if (verbose)
                std::cout << "    " << std::setw(4) << (it+1) << "   " << resNorm
                          << "   " << resNorm/defScale << std::endl;

            result.converged = (resNorm <= reduction * defScale);
            if (std::abs(omega) < 1e-300) break;
        }

        // x is consistent; make sure it stays so for the caller
        copyOwnerToAll_(x);

        result.iterations = static_cast<int>(it);
        result.reduction = resNorm / defScale;
        result.conv_rate = (it > 0) ? std::pow(result.reduction, 1.0/it) : 0.0;
        return result;
    }

    // ---- configuration and local subdomain solver ------------------------------------------

    void readParams_()
    {
        localSolverType_ = getParamFromGroup<std::string>(
            this->paramGroup(), "LinearSolver.Preconditioner.LocalSolver", "mumps");
        twoLevel_ = getParamFromGroup<bool>(
            this->paramGroup(), "LinearSolver.Preconditioner.TwoLevel", false);
        iluOrder_ = getParamFromGroup<int>(
            this->paramGroup(), "LinearSolver.Preconditioner.ILUOrder", 0);
        iluRelaxation_ = getParamFromGroup<double>(
            this->paramGroup(), "LinearSolver.Preconditioner.Relaxation", 1.0);
    }

    //! (Re)build the local block preconditioner for the current local matrix.
    void prepareLocal_(const ScalarMatrix& A)
    {
        if (localSolverType_ == "ilu")
        {
            // inexact local solve: cheap ILU(k), rebuilt each solve (no expensive factorization)
            localPrec_ = std::make_shared<Dune::SeqILU<ScalarMatrix, ScalarBlockVector, ScalarBlockVector>>(
                A, iluOrder_, iluRelaxation_, /*resort=*/false);
        }
        else
        {
            // exact local solve: MUMPS, with the symbolic analysis cached across Newton steps
            if (!mumpsLocal_)
                mumpsLocal_ = std::make_shared<MumpsLocal>(A, this->verbosity() > 1 ? 1 : 0);
            else
                mumpsLocal_->refactorize(A);
            localPrec_ = mumpsLocal_;
        }
    }

    // ---- two-level Schwarz coarse correction -----------------------------------------------

    //! Determine the globally active coarse dofs: one constant per (rank, field, equation
    //! component) whose support (owned dofs of that class on that rank) is non-empty.
    void setupCoarseSpace_()
    {
        std::vector<int> localCount(nLocalCoarse_, 0);
        for (std::size_t i = 0; i < isOwned_.size(); ++i)
            if (isOwned_[i])
                ++localCount[coarseClass_[i]];

        int np = 1, myRank = 0;
        MPI_Comm_size(mpiComm_, &np);
        MPI_Comm_rank(mpiComm_, &myRank);

        std::vector<int> allCount(static_cast<std::size_t>(np) * nLocalCoarse_, 0);
        MPI_Allgather(localCount.data(), static_cast<int>(nLocalCoarse_), MPI_INT,
                      allCount.data(), static_cast<int>(nLocalCoarse_), MPI_INT, mpiComm_);

        coarseRank_.clear();
        coarseLocalClass_.clear();
        myCoarseGlobalIdx_.assign(nLocalCoarse_, -1);
        int g = 0;
        for (int q = 0; q < np; ++q)
            for (std::size_t c = 0; c < nLocalCoarse_; ++c)
                if (allCount[static_cast<std::size_t>(q) * nLocalCoarse_ + c] > 0)
                {
                    coarseRank_.push_back(q);
                    coarseLocalClass_.push_back(static_cast<int>(c));
                    if (q == myRank)
                        myCoarseGlobalIdx_[c] = g;
                    ++g;
                }
        nCoarse_ = g;
    }

    //! Assemble the Galerkin coarse operator A0 = R0 A R0^T (small, replicated on every rank).
    void buildCoarseOperator_(const ScalarMatrix& A)
    {
        const std::size_t nc = static_cast<std::size_t>(nCoarse_);
        std::vector<double> localA0(nc * nc, 0.0);

        ScalarBlockVector phi(isOwned_.size()), w;
        for (int J = 0; J < nCoarse_; ++J)
        {
            // build the consistent coarse basis vector phi_J (indicator of its owned support)
            phi = 0.0;
            if (coarseRank_[J] == rank_)
            {
                const int c = coarseLocalClass_[J];
                for (std::size_t i = 0; i < isOwned_.size(); ++i)
                    if (isOwned_[i] && static_cast<int>(coarseClass_[i]) == c)
                        phi[i] = 1.0;
            }
            copyOwnerToAll_(phi);

            applyA_(A, phi, w); // consistent A phi_J

            // A0[I][J] = phi_I^T (A phi_J) = sum over this rank's owned dofs grouped by class I
            for (std::size_t i = 0; i < isOwned_.size(); ++i)
                if (isOwned_[i])
                {
                    const int I = myCoarseGlobalIdx_[coarseClass_[i]];
                    if (I >= 0)
                        localA0[static_cast<std::size_t>(I) * nc + static_cast<std::size_t>(J)] += w[i][0];
                }
        }

        std::vector<double> fullA0(nc * nc, 0.0);
        MPI_Allreduce(localA0.data(), fullA0.data(), static_cast<int>(nc * nc), MPI_DOUBLE, MPI_SUM, mpiComm_);

        A0_.resize(nc, nc);
        for (std::size_t I = 0; I < nc; ++I)
            for (std::size_t J = 0; J < nc; ++J)
                A0_[I][J] = fullA0[I * nc + J];
    }

    //! Apply the coarse correction zc = R0^T A0^{-1} R0 r for a consistent residual r.
    void applyCoarse_(const ScalarBlockVector& r, ScalarBlockVector& zc) const
    {
        const std::size_t nc = static_cast<std::size_t>(nCoarse_);

        // restrict: rhs0[J] = phi_J^T r  (only the owner rank of J contributes)
        std::vector<double> localRhs(nc, 0.0);
        for (std::size_t i = 0; i < isOwned_.size(); ++i)
            if (isOwned_[i])
            {
                const int I = myCoarseGlobalIdx_[coarseClass_[i]];
                if (I >= 0)
                    localRhs[static_cast<std::size_t>(I)] += r[i][0];
            }
        std::vector<double> rhs0(nc, 0.0);
        MPI_Allreduce(localRhs.data(), rhs0.data(), static_cast<int>(nc), MPI_DOUBLE, MPI_SUM, mpiComm_);

        // solve the small replicated coarse system A0 z0 = rhs0
        Dune::DynamicVector<double> b0(nc), z0(nc, 0.0);
        for (std::size_t J = 0; J < nc; ++J) b0[J] = rhs0[J];
        A0_.solve(z0, b0);

        // prolongate: zc = R0^T z0  (constant z0[J] on the owned support of phi_J), then broadcast
        zc.resize(isOwned_.size());
        zc = 0.0;
        for (std::size_t i = 0; i < isOwned_.size(); ++i)
            if (isOwned_[i])
            {
                const int I = myCoarseGlobalIdx_[coarseClass_[i]];
                if (I >= 0)
                    zc[i] = z0[static_cast<std::size_t>(I)];
            }
        copyOwnerToAll_(zc);
    }

    MPI_Comm mpiComm_;
    int rank_ = 0;
    bool isParallel_;

    std::vector<bool> isOwned_;      // per local scalar dof
    std::vector<SubdomainInfo> subs_;

    // local subdomain solver configuration
    std::string localSolverType_ = "mumps"; // "mumps" (exact) or "ilu" (inexact)
    int iluOrder_ = 0;
    double iluRelaxation_ = 1.0;
    std::shared_ptr<LocalPreconditioner> localPrec_;
    // cached MUMPS local factorization (symbolic analysis reused across Newton steps)
    std::shared_ptr<MumpsLocal> mumpsLocal_;

    // two-level Schwarz coarse correction
    bool twoLevel_ = false;
    std::size_t nLocalCoarse_ = 0;          // coarse classes on this rank (sum of field block sizes)
    std::vector<std::size_t> coarseClass_;  // per local scalar dof: its coarse class
    int nCoarse_ = 0;                       // number of globally active coarse dofs
    std::vector<int> coarseRank_;           // owner rank of each global coarse dof
    std::vector<int> coarseLocalClass_;     // local class of each global coarse dof
    std::vector<int> myCoarseGlobalIdx_;    // this rank's class -> global coarse index (or -1)
    Dune::DynamicMatrix<double> A0_;        // replicated dense coarse operator
};

} // end namespace Dumux

#endif // DUMUX_HAVE_MUMPS

#endif
