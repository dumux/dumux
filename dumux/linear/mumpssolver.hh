// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Direct sparse solver calling MUMPS
 *
 * Supports single-domain (BCRSMatrix/BlockVector) systems and MultiType/multidomain
 * systems (flattened to scalar via MatrixConverter/VectorConverter), sequential and
 * parallel (non-overlapping and overlapping) decompositions.
 *
 * The multidomain support is agnostic to the number of subdomains and allows each
 * subdomain to live on its own grid (e.g. a coupled free-flow / porous-medium problem
 * with momentum, mass and Darcy subdomains on two different grids). The parallel
 * constructor takes a tuple of the subdomain grid geometries (in the same order as the
 * MultiTypeBlockVector subvectors); the sequential path needs no grid information and is
 * N-agnostic by construction (it just flattens and renumbers the global system).
 *
 * This solver backend requires MPI (even for sequential runs, using MPI_COMM_SELF).
 *
 * Note on the RHS: MUMPS also offers a distributed RHS/solution (ICNTL(20)=10/11,
 * ICNTL(21)=1) that would avoid the host gather entirely. In testing it produced
 * incorrect results / crashes when combined with the distributed assembled matrix
 * (ICNTL(18)=3) in this MUMPS 5.7.3 build, so we keep the centralized dense RHS.
 * It is not a bottleneck in practice (O(N) vs the O(nnz) matrix and the factorization).
 */
#ifndef DUMUX_LINEAR_MUMPS_SOLVER_HH
#define DUMUX_LINEAR_MUMPS_SOLVER_HH

#include <string>
#include <vector>
#include <tuple>
#include <numeric>
#include <limits>
#include <algorithm>
#include <functional>
#include <chrono>
#include <iostream>
#include <cstring>

#include <dune/common/exceptions.hh>
#include <dune/common/indices.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <dumux/common/parameters.hh>
#include <dumux/linear/solver.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/matrixconverter.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

#if DUMUX_HAVE_MUMPS

#include <mpi.h>
#include <dmumps_c.h>

namespace Dumux::Detail {

// Provides GridView, DofMapper, dofCodim from LSTraits when canCommunicate=true,
// or inert dummy types when canCommunicate=false (SeqLinearSolverTraits).
template<class LSTraits, bool = LSTraits::canCommunicate>
struct MumpsGridTypes {
    struct GridView {};
    struct DofMapper { std::size_t size() const { return 0; } };
    static constexpr int dofCodim = 0;
};

template<class LSTraits>
struct MumpsGridTypes<LSTraits, true> {
    using GridView = typename LSTraits::GridView;
    using DofMapper = typename LSTraits::DofMapper;
    static constexpr int dofCodim = LSTraits::dofCodim;
};

// Dereference pointer-like objects (e.g. std::shared_ptr<GridGeometry>) and pass through
// references/values unchanged. Used so the multidomain constructor accepts a tuple of either
// (shared) pointers to grid geometries or the grid geometries themselves.
template<class T>
decltype(auto) mumpsDeref(const T& t)
{
    if constexpr (requires (const T& x) { *x; })
        return (*t);
    else
        return (t);
}

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Direct sparse solver calling MUMPS directly via distributed assembled input.
 */
template<class LSTraits, class LATraits>
class DirectSolverMumps : public LinearSolver
{
    using Matrix = typename LATraits::Matrix;
    using XVector = typename LATraits::Vector;
    using BVector = typename LATraits::Vector;

    using Scalar = double;
    // Global scalar-DOF index. MUMPS itself uses 32-bit MUMPS_INT for the matrix
    // entries; we keep a wider type internally and cast (with a check) at the boundary.
    using GlobalIndex = long long;

    using GridTypes_ = Detail::MumpsGridTypes<LSTraits>;
    using GridView = typename GridTypes_::GridView;
    using DofMapper = typename GridTypes_::DofMapper;
    static constexpr int dofCodim = GridTypes_::dofCodim;

    // True when XVector is a MultiTypeBlockVector (multidomain assembler).
    static constexpr bool isMultiType = isMultiTypeBlockVector<XVector>::value;
    // Scalar DOFs per block DOF. MultiType is always converted to scalar (blockSize=1).
    static constexpr std::size_t blockSize = []() constexpr -> std::size_t {
        if constexpr (isMultiTypeBlockVector<XVector>::value) return 1;
        else return XVector::block_type::size();
    }();

    // Flattened scalar block vector type produced by VectorConverter (1 scalar per block).
    // The per-subdomain RHS/ghost closures operate on slices of a vector of this type.
    using ScalarBlockVector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;

    // Per-subdomain parallel bookkeeping (MultiType only). One entry per subdomain, in the
    // same order as the MultiTypeBlockVector subvectors / the flattened scalar layout.
    struct SubdomainInfo
    {
        std::size_t scalarBS = 1;        // scalar DOFs per block DOF in this subdomain
        std::size_t nBlockLocal = 0;     // local block DOFs of this subdomain on this rank
        std::size_t scalarOffset = 0;    // start index of this subdomain in the flattened scalar vector
        bool rowsComplete = true;        // owner border rows already complete after local assembly?
        // Sum partial border-DOF RHS contributions onto the owner (non-overlapping incomplete rows).
        // Empty when not needed (overlapping or complete rows). Operates on its scalar slice.
        std::function<void(ScalarBlockVector&)> rhsPreprocess;
        // Broadcast the owner's solution value to all copies (zero non-owned beforehand).
        // Operates on its scalar slice.
        std::function<void(ScalarBlockVector&)> ghostBroadcast;
    };

public:
    struct SolverResult : public Dune::InverseOperatorResult
    {
        SolverResult() = default;
        SolverResult(const Dune::InverseOperatorResult& o) : InverseOperatorResult(o) {}
        operator bool() const { return this->converged; }
    };

    // Constructor for sequential (SeqLinearSolverTraits) solvers with no grid.
    // DOF maps are built lazily on the first solve() call (single-type and MultiType).
    explicit DirectSolverMumps(const std::string& paramGroup = "")
    requires (!LSTraits::canCommunicate)
    : LinearSolver(paramGroup)
    , isParallel_(false)
    {
        mpiComm_ = MPI_COMM_SELF;
        initMumps_();
        readSolverParams_();
    }

    // Constructor for parallel/sequential single-domain solvers with a grid.
    template<class GridGeometry>
    DirectSolverMumps(const GridGeometry& gridGeometry,
                      const GridView& gridView,
                      const DofMapper& dofMapper,
                      const std::string& paramGroup = "")
    requires (LSTraits::canCommunicate)
    : LinearSolver(paramGroup)
    , gridView_(gridView)
    , mapper_(&dofMapper)
    , isParallel_(gridView.comm().size() > 1)
    {
        mpiComm_ = static_cast<MPI_Comm>(gridView.comm());
        initMumps_();
        readSolverParams_();

        buildGlobalDofIndices_(gridView, dofMapper);
    }

    /*!
     * Constructor for parallel MultiType (multidomain) problems with an arbitrary number of
     * subdomains, each possibly on its own grid. \a ggTuple is a tuple of the subdomain grid
     * geometries (or (shared) pointers to them), in the same order as the MultiTypeBlockVector
     * subvectors. Builds the combined parallel DOF maps and the per-subdomain RHS/ghost closures.
     */
    template<class GGTuple>
    explicit DirectSolverMumps(const GGTuple& ggTuple, const std::string& paramGroup = "")
    requires (isMultiType && requires { std::tuple_size<GGTuple>::value; })
    : LinearSolver(paramGroup)
    , isParallel_(Detail::mumpsDeref(std::get<0>(ggTuple)).gridView().comm().size() > 1)
    {
        constexpr std::size_t numSubDomains = std::tuple_size_v<GGTuple>;
        static_assert(numSubDomains == std::tuple_size_v<XVector>,
                      "Number of grid geometries must match the number of MultiTypeBlockVector subdomains");

        mpiComm_ = static_cast<MPI_Comm>(Detail::mumpsDeref(std::get<0>(ggTuple)).gridView().comm());
        initMumps_();
        readSolverParams_();

        subs_.resize(numSubDomains);

        // Running offsets while we append the subdomains' scalar DOFs in order.
        GlobalIndex scalarGlobalOffset = 0;  // global scalar index where the current subdomain starts
        std::size_t scalarLocalOffset = 0;   // local  scalar index where the current subdomain starts

        using Dune::Hybrid::forEach;
        forEach(std::make_index_sequence<numSubDomains>{}, [&](auto I)
        {
            constexpr std::size_t i = decltype(I)::value;
            const auto& gg = Detail::mumpsDeref(std::get<i>(ggTuple));
            using SubGG = std::decay_t<decltype(gg)>;
            using SubTraits = LinearSolverTraits<SubGG>;
            using SubBlockVec = std::tuple_element_t<i, XVector>;
            constexpr std::size_t scalarBS = SubBlockVec::block_type::size();
            using GV = typename SubGG::GridView;
            constexpr bool canComm = subCanComm_<SubTraits, GV>();

            // Parallel runs require every subdomain's DOF mapper to support the generic grid
            // communication; bail out with a clear message otherwise (e.g. FCStaggered momentum).
            if (isParallel_ && !canComm)
                DUNE_THROW(Dune::NotImplemented,
                    "DirectSolverMumps: parallel multidomain solve is not supported for subdomain "
                    << i << " — its discretization's DOF mapper (e.g. FaceCenteredStaggered) cannot "
                    "drive the generic grid communication. This combination works sequentially; for "
                    "parallel runs use a momentum discretization with a full mapper (CVFE/FCDiamond).");

            // Block-level global IDs for this subdomain (numbered within the subdomain, [0, nGlobalBlock)).
            std::vector<GlobalIndex> gids;
            std::vector<bool> owned, ghost;
            const GlobalIndex nGlobalBlock = buildSubdomainGIDs_<SubTraits, canComm>(
                gg.gridView(), subDofMapper_<SubTraits>(gg), 0, gids, owned, ghost);
            const std::size_t nBlockLocal = gids.size();

            // Are the owner's border rows already complete after local (interior-element) assembly?
            // Reproduces the original momentum/mass behaviour generically:
            //  - overlapping decompositions: complete (the overlap supplies neighbour contributions);
            //  - non-overlapping: incomplete for multi-codim schemes (e.g. CVFE/PQ1Bubble momentum)
            //    and for vertex-based schemes (Box), complete otherwise (e.g. face-based FCDiamond mass).
            const bool nonOverlapping = SubTraits::isNonOverlapping(gg.gridView());
            constexpr bool hasMultiCodim = requires { SubTraits::dofCodims; };
            constexpr bool vertexBased = (SubTraits::dofCodim == static_cast<int>(GV::dimension));
            const bool rowsComplete = nonOverlapping ? !(hasMultiCodim || vertexBased) : true;

            // Append this subdomain's scalar DOFs to the global bookkeeping (in subvector order,
            // matching MatrixConverter/VectorConverter's flattened scalar layout).
            for (std::size_t b = 0; b < nBlockLocal; ++b)
                for (std::size_t k = 0; k < scalarBS; ++k)
                {
                    localToGlobal_.push_back(scalarGlobalOffset + scalarBS * gids[b] + static_cast<GlobalIndex>(k));
                    isOwned_.push_back(owned[b]);
                    isGhost_.push_back(ghost[b]);
                    rowComplete_.push_back(rowsComplete);
                }

            subs_[i].scalarBS = scalarBS;
            subs_[i].nBlockLocal = nBlockLocal;
            subs_[i].scalarOffset = scalarLocalOffset;
            subs_[i].rowsComplete = rowsComplete;
            if (isParallel_)
                if constexpr (canComm)
                    setupSubdomainClosures_<i>(gg);

            scalarLocalOffset += scalarBS * nBlockLocal;
            scalarGlobalOffset += static_cast<GlobalIndex>(scalarBS) * nGlobalBlock;
        });

        numBlocksGlobal_ = static_cast<std::size_t>(scalarGlobalOffset);
    }

    ~DirectSolverMumps()
    {
        if (initialized_)
        {
            id_.job = -2; // terminate the MUMPS instance
            dmumps_c(&id_);
        }
    }

    // Non-copyable: owns a MUMPS instance and raw pointers into our buffers.
    DirectSolverMumps(const DirectSolverMumps&) = delete;
    DirectSolverMumps& operator=(const DirectSolverMumps&) = delete;

    SolverResult solve(const Matrix& A, XVector& x, const BVector& b)
    {
        if constexpr (isMultiType)
            return solveMultiType_(A, x, b);
        else
            return solveSingleType_(A, x, b);
    }

    Scalar norm(const XVector& x) const
    {
        if constexpr (isMultiType)
        {
            if (!isParallel_ || localToGlobal_.empty())
            {
                auto y = VectorConverter<XVector>::multiTypeToBlockVector(x);
                return y.two_norm();
            }
            // Parallel: sum border-DOF partial contributions before the owned-DOF global norm.
            auto v = VectorConverter<XVector>::multiTypeToBlockVector(x);
            for (const auto& sub : subs_)
                if (sub.rhsPreprocess)
                    sub.rhsPreprocess(v);

            double localSumSq = 0.0;
            for (std::size_t r = 0; r < v.size(); ++r)
                if (isOwned_[r])
                    localSumSq += v[r][0] * v[r][0];

            double globalSumSq = 0.0;
            MPI_Allreduce(&localSumSq, &globalSumSq, 1, MPI_DOUBLE, MPI_SUM, mpiComm_);
            using std::sqrt;
            return sqrt(globalSumSq);
        }
        else
        {
            if (!isParallel_ || localToGlobal_.empty())
                return x.two_norm();

            // Parallel single-type: true global residual norm over owned DOFs.
            if constexpr (LSTraits::canCommunicate)
            {
                if (LSTraits::isNonOverlapping(gridView_))
                {
                    auto xcopy = x;
                    ParallelVectorHelper<GridView, DofMapper, dofCodim> vecHelper(gridView_, *mapper_);
                    if constexpr (requires { LSTraits::dofCodims; })
                        vecHelper.makeNonOverlappingConsistent(xcopy, LSTraits::dofCodims);
                    else
                        vecHelper.makeNonOverlappingConsistent(xcopy);
                    return ownedGlobalNorm_(xcopy);
                }
            }
            return ownedGlobalNorm_(x);
        }
    }

    std::string name() const { return "DuMux direct MUMPS solver"; }

private:
    // Accessor for MUMPS control parameters using the manual's 1-based ICNTL(i) numbering.
    int& icntl_(int i) { return id_.icntl[i - 1]; }
    const int& infog_(int i) const { return id_.infog[i - 1]; }

    /*!
     * Initialize the MUMPS instance (JOB=-1) and set our default control parameters.
     * par=1: the host also participates in the factorization. sym=0: unsymmetric.
     * ICNTL(18)=3: distributed assembled input.
     */
    void initMumps_()
    {
        // Zero the whole control struct so distributed-RHS/solution fields MUMPS may read
        // (nz_rhs, nloc_rhs, pointers, ...) start well-defined rather than holding stack garbage.
        std::memset(&id_, 0, sizeof(id_));
        id_.comm_fortran = static_cast<int>(MPI_Comm_c2f(mpiComm_));
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
        icntl_(18) = 3; // distributed assembled input: each rank supplies its local triplets
        icntl_(20) = 0; // dense, centralized RHS on the host
        icntl_(21) = 0; // centralized solution on the host
    }

    /*!
     * Read optional MUMPS tuning parameters. LinearSolver.MumpsOrdering sets ICNTL(7),
     * the fill-reducing ordering used during the (reused) symbolic analysis. Values:
     * 0=AMD, 2=AMF, 4=PORD, 6=QAMD, 7=automatic. 3=SCOTCH/5=METIS require those libraries.
     */
    void readSolverParams_()
    {
        const int ordering = getParamFromGroup<int>(this->paramGroup(), "LinearSolver.MumpsOrdering", -1);
        if (ordering >= 0)
            icntl_(7) = ordering;

        // ICNTL(28): symbolic analysis mode. 0=automatic, 1=sequential (default), 2=parallel.
        // Parallel analysis requires MUMPS built with a parallel ordering tool (PT-SCOTCH or
        // ParMETIS); without one MUMPS silently falls back to sequential analysis. ICNTL(29)
        // selects that tool (0=auto, 1=PT-SCOTCH, 2=ParMETIS).
        const int parAnalysis = getParamFromGroup<int>(this->paramGroup(), "LinearSolver.MumpsParallelAnalysis", -1);
        if (parAnalysis >= 0)
            icntl_(28) = parAnalysis;
        const int parOrderingTool = getParamFromGroup<int>(this->paramGroup(), "LinearSolver.MumpsParallelOrderingTool", -1);
        if (parOrderingTool >= 0)
            icntl_(29) = parOrderingTool;
    }

    /*!
     * Compute localToGlobal_[i] = global block index for Dune local block DOF i,
     * plus isOwned_[i] and isGhost_[i] (single-domain case).
     */
    void buildGlobalDofIndices_(const GridView& gridView, const DofMapper& mapper)
    {
        const std::size_t nLocal = mapper.size();
        localToGlobal_.assign(nLocal, std::numeric_limits<GlobalIndex>::max());
        isOwned_.assign(nLocal, false);
        isGhost_.assign(nLocal, false);

        if constexpr (LSTraits::canCommunicate)
        {
            if (isParallel_)
            {
                ParallelISTLHelper<LSTraits> pHelper(gridView, mapper);
                for (std::size_t i = 0; i < nLocal; ++i)
                {
                    isOwned_[i] = pHelper.isOwned(i);
                    isGhost_[i] = pHelper.isGhost(i);
                }

                GlobalIndex numOwned = static_cast<GlobalIndex>(
                    std::count(isOwned_.begin(), isOwned_.end(), true));
                const int commSize = gridView.comm().size();
                const int myRank = gridView.comm().rank();
                std::vector<GlobalIndex> counts(commSize);
                gridView.comm().allgather(&numOwned, 1, counts.data());

                GlobalIndex offset = 0;
                for (int p = 0; p < myRank; ++p)
                    offset += counts[p];
                numBlocksGlobal_ = 0;
                for (auto c : counts) numBlocksGlobal_ += static_cast<std::size_t>(c);

                for (std::size_t i = 0; i < nLocal; ++i)
                    if (isOwned_[i])
                        localToGlobal_[i] = offset++;

                // Propagate global IDs from owners to ghost DOFs via min communication.
                if constexpr (requires { LSTraits::dofCodims; }) {
                    MultiCodimVectorCommDataHandleMin<DofMapper, std::vector<GlobalIndex>, GridView::dimension, GlobalIndex>
                        handle(mapper, localToGlobal_, LSTraits::dofCodims);
                    gridView.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
                } else {
                    VectorCommDataHandleMin<DofMapper, std::vector<GlobalIndex>, dofCodim, GlobalIndex>
                        handle(mapper, localToGlobal_);
                    gridView.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
                }
                return;
            }
        }

        // Sequential: all DOFs owned with consecutive IDs.
        numBlocksGlobal_ = nLocal;
        std::fill(isOwned_.begin(), isOwned_.end(), true);
        for (std::size_t i = 0; i < nLocal; ++i)
            localToGlobal_[i] = static_cast<GlobalIndex>(i);
    }

    // Obtain a subdomain's DOF mapper discretization-generically: most grid geometries expose
    // dofMapper(); FCStaggered only provides one through its linear-solver traits (a lightweight
    // by-value wrapper around the grid view).
    template<class SubTraits, class SubGG>
    static decltype(auto) subDofMapper_(const SubGG& gg)
    {
        if constexpr (requires { gg.dofMapper(); })
            return (gg.dofMapper());
        else
            return SubTraits::dofMapper(gg);
    }

    // True if a subdomain's DOF mapper can drive the generic codim-based grid communication
    // (it must expose indices(entity)); the minimal FCStaggered mapper wrapper cannot, so
    // FCStaggered subdomains are only supported sequentially.
    template<class SubTraits, class GV>
    static constexpr bool subCanComm_()
    {
        using Mapper = typename SubTraits::DofMapper;
        using Element = typename GV::template Codim<0>::Entity;
        return SubTraits::canCommunicate
            && requires (const Mapper& m, const Element& e) { m.indices(e); };
    }

    // Build global block indices for one subdomain (MultiType). Returns the subdomain's
    // total number of global block DOFs. gids[i] = global block ID (offset to be globally unique).
    // CanComm gates the parallel branch so it is not instantiated for mappers that cannot
    // communicate (e.g. FCStaggered); such subdomains fall back to sequential numbering.
    template<class SubLSTraits, bool CanComm, class GV, class Mapper>
    GlobalIndex buildSubdomainGIDs_(const GV& gridView, const Mapper& mapper,
                                    GlobalIndex globalOffset,
                                    std::vector<GlobalIndex>& gids,
                                    std::vector<bool>& owned,
                                    std::vector<bool>& ghost)
    {
        const std::size_t nLocal = mapper.size();
        gids.assign(nLocal, std::numeric_limits<GlobalIndex>::max());
        owned.assign(nLocal, false);
        ghost.assign(nLocal, false);

        if constexpr (CanComm)
        {
            if (isParallel_)
            {
                ParallelISTLHelper<SubLSTraits> pHelper(gridView, mapper);
                for (std::size_t i = 0; i < nLocal; ++i)
                {
                    owned[i] = pHelper.isOwned(i);
                    ghost[i] = pHelper.isGhost(i);
                }

                GlobalIndex numOwned = static_cast<GlobalIndex>(
                    std::count(owned.begin(), owned.end(), true));
                const int commSize = gridView.comm().size();
                const int myRank = gridView.comm().rank();
                std::vector<GlobalIndex> counts(commSize);
                gridView.comm().allgather(&numOwned, 1, counts.data());

                GlobalIndex localOffset = globalOffset;
                for (int p = 0; p < myRank; ++p)
                    localOffset += counts[p];

                GlobalIndex nGlobal = 0;
                for (auto c : counts) nGlobal += c;

                for (std::size_t i = 0; i < nLocal; ++i)
                    if (owned[i])
                        gids[i] = localOffset++;

                if constexpr (requires { SubLSTraits::dofCodims; }) {
                    MultiCodimVectorCommDataHandleMin<Mapper, std::vector<GlobalIndex>, GV::dimension, GlobalIndex>
                        handle(mapper, gids, SubLSTraits::dofCodims);
                    gridView.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
                } else {
                    constexpr int codim = SubLSTraits::dofCodim;
                    VectorCommDataHandleMin<Mapper, std::vector<GlobalIndex>, codim, GlobalIndex>
                        handle(mapper, gids);
                    gridView.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
                }
                return nGlobal;
            }
        }

        // Sequential subdomain: all DOFs owned with consecutive IDs.
        std::fill(owned.begin(), owned.end(), true);
        for (std::size_t i = 0; i < nLocal; ++i)
            gids[i] = globalOffset + static_cast<GlobalIndex>(i);
        return static_cast<GlobalIndex>(nLocal);
    }

    /*!
     * Build the RHS-preprocessing and ghost-broadcast closures for subdomain \a i (parallel only).
     * The closures operate on the subdomain's contiguous slice [scalarOffset, scalarOffset +
     * scalarBS*nBlockLocal) of the flattened scalar vector: they unpack it into a typed subdomain
     * block vector, run the Dune grid communication, and pack the result back. Type-erasing into
     * std::function lets us hold heterogeneous subdomains (different grids, block sizes, codims) in
     * one container.
     *
     * No Tpetra-style extra-count communication is needed: MUMPS sums the duplicate border entries
     * supplied by the different ranks during the distributed assembly.
     */
    template<std::size_t i, class SubGG>
    void setupSubdomainClosures_(const SubGG& gg)
    {
        using SubTraits = LinearSolverTraits<SubGG>;
        using SubBlockVec = std::tuple_element_t<i, XVector>;
        using GV = typename SubGG::GridView;
        using Mapper = typename SubTraits::DofMapper;
        constexpr std::size_t scalarBS = SubBlockVec::block_type::size();

        const GV gv = gg.gridView();
        // Store our own copy of the DOF mapper so the captured closures own a stable instance
        // regardless of whether the grid geometry hands it out by reference or by value.
        auto mapper = std::make_shared<Mapper>(subDofMapper_<SubTraits>(gg));
        const std::size_t offset = subs_[i].scalarOffset;
        const std::size_t nBlock = subs_[i].nBlockLocal;

        // Pack/unpack helpers between the scalar slice and a typed subdomain block vector.
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

        // Ghost broadcast: sum each DOF's copies to the owner's value (callers zero non-owned first).
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

        // RHS preprocessing: only needed for non-overlapping decompositions with incomplete owner
        // rows (sum partial border contributions onto the owner before the centralized gather).
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

    // Initialize DOF maps for sequential use from a known scalar DOF count.
    // Called lazily on first solve for the no-grid constructor (single-type and MultiType).
    void initSeqMaps_(std::size_t nDofs)
    {
        numBlocksGlobal_ = nDofs;
        localToGlobal_.resize(nDofs);
        isOwned_.assign(nDofs, true);
        isGhost_.assign(nDofs, false);
        rowComplete_.assign(nDofs, true);
        for (std::size_t i = 0; i < nDofs; ++i)
            localToGlobal_[i] = static_cast<GlobalIndex>(i);
    }

    /*!
     * Mark the block rows this rank contributes to the distributed assembly.
     *  - MultiType: per scalar row, isOwned_ if the owner row is already complete, else !isGhost_
     *    (so every rank supplies its partial border row and MUMPS sums them).
     *  - single-domain non-overlapping: every non-ghost row (owner border rows incomplete).
     *  - single-domain overlapping: only owned rows (owner rows already complete).
     *  - sequential: all owned rows.
     */
    void buildIsSrcRow_()
    {
        const std::size_t n = isOwned_.size();
        isSrcRow_.assign(n, false);
        if constexpr (isMultiType)
        {
            for (std::size_t r = 0; r < n; ++r)
                isSrcRow_[r] = rowComplete_[r] ? isOwned_[r] : !isGhost_[r];
            return;
        }
        else if constexpr (LSTraits::canCommunicate)
        {
            if (isParallel_)
            {
                const bool nonOverlapping = LSTraits::isNonOverlapping(gridView_);
                for (std::size_t r = 0; r < n; ++r)
                    isSrcRow_[r] = nonOverlapping ? !isGhost_[r] : isOwned_[r];
                return;
            }
        }
        isSrcRow_ = isOwned_;
    }

    /*!
     * Build the distributed triplet structure (irn_loc/jcn_loc) once, in a fixed iteration
     * order over the contributing rows. a_loc is (re)filled in the same order on every solve,
     * so the symbolic analysis stays valid. Indices are global and 1-based (Fortran).
     */
    template<class AMatrix>
    void assembleTriplets_(const AMatrix& A, bool firstCall)
    {
        using BlockType = typename AMatrix::block_type;
        constexpr int BS_I = BlockType::rows;
        constexpr int BS_J = BlockType::cols;

        if (firstCall)
        {
            std::size_t nnzLocal = 0;
            for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
                if (isSrcRow_[rowIt.index()])
                    nnzLocal += rowIt->size() * BS_I * BS_J;
            irnLoc_.resize(nnzLocal);
            jcnLoc_.resize(nnzLocal);
            aLoc_.resize(nnzLocal);
        }

        std::size_t k = 0;
        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            const std::size_t r = rowIt.index();
            if (!isSrcRow_[r]) continue;
            const GlobalIndex gBlockRow = localToGlobal_[r];

            for (int ki = 0; ki < BS_I; ++ki)
            {
                const GlobalIndex gRow = gBlockRow * BS_I + ki;
                for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                {
                    const GlobalIndex gBlockCol = localToGlobal_[colIt.index()];
                    for (int kj = 0; kj < BS_J; ++kj)
                    {
                        if (firstCall)
                        {
                            irnLoc_[k] = toMumpsInt_(gRow + 1);                  // 1-based
                            jcnLoc_[k] = toMumpsInt_(gBlockCol * BS_J + kj + 1); // 1-based
                        }
                        aLoc_[k] = (*colIt)[ki][kj];
                        ++k;
                    }
                }
            }
        }
    }

    /*!
     * One-time gather metadata for the centralized (host) dense RHS/solution. Each rank's owned
     * scalar DOFs may map to a non-contiguous set of global indices (e.g. MultiType: momentum
     * DOFs are numbered globally before mass DOFs), so we record the global index of every owned
     * scalar value and use it to permute between gather order and natural (MUMPS) order.
     */
    void setupRhsGather_()
    {
        // Owned scalar global indices in the (ascending-local) packing order.
        std::vector<int> myIdx;
        for (std::size_t i = 0; i < localToGlobal_.size(); ++i)
            if (isOwned_[i])
                for (std::size_t k = 0; k < blockSize; ++k)
                    myIdx.push_back(static_cast<int>(localToGlobal_[i] * blockSize + k));
        nOwnedScalar_ = static_cast<int>(myIdx.size());

        int commSize = 1, myRank = 0;
        MPI_Comm_size(mpiComm_, &commSize);
        MPI_Comm_rank(mpiComm_, &myRank);
        isHost_ = (myRank == 0);

        recvCounts_.assign(isHost_ ? commSize : 0, 0);
        MPI_Gather(&nOwnedScalar_, 1, MPI_INT,
                   isHost_ ? recvCounts_.data() : nullptr, 1, MPI_INT, 0, mpiComm_);

        const std::size_t nScalar = static_cast<std::size_t>(numBlocksGlobal_ * blockSize);
        if (isHost_)
        {
            displs_.assign(commSize, 0);
            std::partial_sum(recvCounts_.begin(), recvCounts_.end() - 1, displs_.begin() + 1);
            allGlobalIdx_.assign(nScalar, 0);
            rhsHost_.assign(nScalar, 0.0);   // natural (MUMPS) order, length n
            tmpHost_.assign(nScalar, 0.0);   // gather order
        }
        // Gather the owned global indices (constant across solves).
        MPI_Gatherv(myIdx.data(), nOwnedScalar_, MPI_INT,
                    isHost_ ? allGlobalIdx_.data() : nullptr,
                    isHost_ ? recvCounts_.data() : nullptr,
                    isHost_ ? displs_.data() : nullptr,
                    MPI_INT, 0, mpiComm_);
        sendBuf_.resize(static_cast<std::size_t>(nOwnedScalar_));
    }

    /*!
     * Gather the owned RHS entries onto the host as a dense vector in natural global order, run
     * the MUMPS solve (JOB=3), and scatter the solution back to the owned DOFs of x.
     */
    template<class BVec, class XVec>
    void solveCentralized_(const BVec& b, XVec& x)
    {
        // Pack owned scalar RHS values in the same order as setupRhsGather_'s indices.
        std::size_t s = 0;
        for (std::size_t i = 0; i < b.size(); ++i)
            if (isOwned_[i])
                for (std::size_t k = 0; k < blockSize; ++k)
                    sendBuf_[s++] = b[i][k];

        MPI_Gatherv(sendBuf_.data(), nOwnedScalar_, MPI_DOUBLE,
                    isHost_ ? tmpHost_.data() : nullptr,
                    isHost_ ? recvCounts_.data() : nullptr,
                    isHost_ ? displs_.data() : nullptr,
                    MPI_DOUBLE, 0, mpiComm_);

        if (isHost_)
        {
            for (std::size_t p = 0; p < tmpHost_.size(); ++p)
                rhsHost_[static_cast<std::size_t>(allGlobalIdx_[p])] = tmpHost_[p];
            id_.nrhs = 1;
            id_.lrhs = id_.n;
            id_.rhs = rhsHost_.data();
        }
        id_.job = 3; // solve
        dmumps_c(&id_);
        checkError_("solve");

        // Solution is centralized on the host (natural order); permute back to gather order.
        if (isHost_)
            for (std::size_t p = 0; p < tmpHost_.size(); ++p)
                tmpHost_[p] = rhsHost_[static_cast<std::size_t>(allGlobalIdx_[p])];

        MPI_Scatterv(isHost_ ? tmpHost_.data() : nullptr,
                     isHost_ ? recvCounts_.data() : nullptr,
                     isHost_ ? displs_.data() : nullptr,
                     MPI_DOUBLE, sendBuf_.data(), nOwnedScalar_, MPI_DOUBLE, 0, mpiComm_);

        s = 0;
        for (std::size_t i = 0; i < x.size(); ++i)
            if (isOwned_[i])
                for (std::size_t k = 0; k < blockSize; ++k)
                    x[i][k] = sendBuf_[s++];
    }

    // Core MUMPS pipeline shared by single-type and MultiType: analysis (once, reused) +
    // numeric factorization + centralized solve. A is the scalar matrix.
    template<class AScalar, class BScalar, class XScalar>
    void solveCore_(const AScalar& A, const BScalar& b, XScalar& x)
    {
        using Clock = std::chrono::steady_clock;
        const bool time = this->verbosity() > 0;
        auto ms = [](Clock::time_point a, Clock::time_point z) {
            return std::chrono::duration<double, std::milli>(z - a).count();
        };

        double tAssemble = 0.0, tAnalysis = 0.0;
        auto t0 = Clock::now();
        if (!analyzed_)
        {
            buildIsSrcRow_();
            setupRhsGather_();

            assembleTriplets_(A, /*firstCall=*/true);
            id_.n = toMumpsInt_(static_cast<GlobalIndex>(numBlocksGlobal_ * blockSize));
            id_.nz_loc = toMumpsInt_(static_cast<GlobalIndex>(aLoc_.size()));
            id_.irn_loc = irnLoc_.data();
            id_.jcn_loc = jcnLoc_.data();
            id_.a_loc = aLoc_.data();
            auto t1 = Clock::now(); tAssemble = ms(t0, t1);

            id_.job = 1; // symbolic analysis (reused across solves)
            dmumps_c(&id_);
            checkError_("symbolic analysis");
            analyzed_ = true;
            tAnalysis = ms(t1, Clock::now());
        }
        else
        {
            assembleTriplets_(A, /*firstCall=*/false);
            id_.a_loc = aLoc_.data();
            tAssemble = ms(t0, Clock::now());
        }

        auto t2 = Clock::now();
        id_.job = 2; // numeric factorization
        dmumps_c(&id_);
        checkError_("numeric factorization");
        const double tNumeric = ms(t2, Clock::now());

        auto t3 = Clock::now();
        solveCentralized_(b, x);
        const double tSolve = ms(t3, Clock::now());

        if (time && isHost_)
            std::cout << "[MUMPS] assemble " << tAssemble << " ms"
                      << (tAnalysis > 0.0 ? (", analysis " + std::to_string(tAnalysis) + " ms") : std::string{})
                      << ", numeric factorization " << tNumeric << " ms"
                      << ", solve+gather " << tSolve << " ms" << std::endl;
    }

    SolverResult solveSingleType_(const Matrix& A, XVector& x, const BVector& b)
    {
        auto bcopy = b;

        // Non-overlapping parallel: accumulate partial border RHS onto owners before the gather.
        if constexpr (LSTraits::canCommunicate)
        {
            if (isParallel_ && LSTraits::isNonOverlapping(gridView_))
            {
                ParallelVectorHelper<GridView, DofMapper, dofCodim> vecHelper(gridView_, *mapper_);
                if constexpr (requires { LSTraits::dofCodims; })
                    vecHelper.makeNonOverlappingConsistent(bcopy, LSTraits::dofCodims);
                else
                    vecHelper.makeNonOverlappingConsistent(bcopy);
            }
        }

        if (localToGlobal_.empty())
            initSeqMaps_(A.N());

        solveCore_(A, bcopy, x);

        // Distribute the solution: zero non-owned, then sum-communicate to the owner's value.
        if constexpr (LSTraits::canCommunicate)
        {
            if (isParallel_)
            {
                for (std::size_t i = 0; i < x.size(); ++i)
                    if (!isOwned_[i])
                        x[i] = 0;
                if constexpr (requires { LSTraits::dofCodims; }) {
                    MultiCodimVectorCommDataHandleSum<DofMapper, XVector, GridView::dimension> handle(*mapper_, x, LSTraits::dofCodims);
                    gridView_.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
                } else {
                    VectorCommDataHandleSum<DofMapper, XVector, dofCodim> handle(*mapper_, x);
                    gridView_.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
                }
            }
        }

        SolverResult result;
        result.converged = true;
        return result;
    }

    SolverResult solveMultiType_(const Matrix& A, XVector& x, const BVector& b)
    {
        // Flatten MultiTypeBlockMatrix/Vector to scalar BCRSMatrix/BlockVector.
        auto Ascalar = MatrixConverter<Matrix>::multiTypeToBCRSMatrix(A);
        auto bscalar = VectorConverter<BVector>::multiTypeToBlockVector(b);

        if (localToGlobal_.empty())
            initSeqMaps_(Ascalar.N());

        // RHS preprocessing per subdomain: sum partial border-DOF contributions onto the owner.
        if (isParallel_)
            for (const auto& sub : subs_)
                if (sub.rhsPreprocess)
                    sub.rhsPreprocess(bscalar);

        auto xscalar = bscalar; // same layout; owned entries overwritten by the solve
        solveCore_(Ascalar, bscalar, xscalar);

        // Ghost broadcast per subdomain: zero non-owned DOFs, sum-communicate to the owner's value.
        if (isParallel_)
        {
            for (std::size_t r = 0; r < xscalar.size(); ++r)
                if (!isOwned_[r]) xscalar[r] = 0.0;
            for (const auto& sub : subs_)
                if (sub.ghostBroadcast)
                    sub.ghostBroadcast(xscalar);
        }

        VectorConverter<XVector>::retrieveValues(x, xscalar);

        SolverResult result;
        result.converged = true;
        return result;
    }

    Scalar ownedGlobalNorm_(const XVector& x) const
    {
        double localSumSq = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i)
            if (isOwned_[i])
                for (std::size_t k = 0; k < blockSize; ++k)
                    localSumSq += x[i][k] * x[i][k];
        double globalSumSq = 0.0;
        MPI_Allreduce(&localSumSq, &globalSumSq, 1, MPI_DOUBLE, MPI_SUM, mpiComm_);
        using std::sqrt;
        return sqrt(globalSumSq);
    }

    static MUMPS_INT toMumpsInt_(GlobalIndex v)
    {
        if (v > static_cast<GlobalIndex>(std::numeric_limits<MUMPS_INT>::max()))
            DUNE_THROW(Dune::RangeError, "Index " << v << " exceeds MUMPS_INT range; "
                       "rebuild MUMPS with 64-bit integers (-DINTSIZE64)");
        return static_cast<MUMPS_INT>(v);
    }

    void checkError_(const char* phase) const
    {
        if (infog_(1) < 0)
            DUNE_THROW(Dune::Exception, "MUMPS " << phase << " failed: INFOG(1)="
                       << infog_(1) << ", INFOG(2)=" << infog_(2));
    }

    GridView gridView_;
    const DofMapper* mapper_ = nullptr;
    bool isParallel_;

    MPI_Comm mpiComm_;
    DMUMPS_STRUC_C id_;
    bool initialized_ = false;
    bool analyzed_ = false;

    // Global DOF bookkeeping (block level for single-type, scalar level for MultiType).
    std::size_t numBlocksGlobal_ = 0;
    std::vector<GlobalIndex> localToGlobal_;
    std::vector<bool> isOwned_;
    std::vector<bool> isGhost_;
    std::vector<bool> isSrcRow_; // block rows this rank contributes to the distributed assembly
    std::vector<bool> rowComplete_; // per (scalar) row: owner row already complete? (MultiType)

    // Distributed assembled input (ICNTL(18)=3). irn_loc/jcn_loc are fixed after the first
    // assembly (so the symbolic analysis is reused); a_loc is refilled every solve.
    std::vector<MUMPS_INT> irnLoc_;
    std::vector<MUMPS_INT> jcnLoc_;
    std::vector<double> aLoc_;

    // Centralized dense RHS/solution gather (host = rank 0). The assembled RHS is gathered to the
    // host, solved, and the solution scattered back. (MUMPS's distributed RHS feature, ICNTL(20)=10,
    // proved unreliable with this build's ICNTL(18)=3 distributed matrix; see the class docs.)
    bool isHost_ = true;
    int nOwnedScalar_ = 0;
    std::vector<int> recvCounts_;     // host only
    std::vector<int> displs_;         // host only
    std::vector<int> allGlobalIdx_;   // host only, gather-order -> natural global index
    std::vector<double> rhsHost_;     // host only, natural (MUMPS) order, length n
    std::vector<double> tmpHost_;     // host only, gather order
    std::vector<double> sendBuf_;     // per-rank owned scalar values / solution

    // MultiType (multidomain) bookkeeping: one entry per subdomain, in subvector order.
    std::vector<SubdomainInfo> subs_;
};

} // end namespace Dumux

#endif // DUMUX_HAVE_MUMPS

#endif
