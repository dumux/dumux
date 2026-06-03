// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Linear solvers from Trilinos
 * This solver backend requires MPI (even for sequential runs)
 */
#ifndef DUMUX_LINEAR_TRILINOS_SOLVERS_HH
#define DUMUX_LINEAR_TRILINOS_SOLVERS_HH

#include <string>
#include <numeric>
#include <limits>
#include <algorithm>
#include <functional>

#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/linear/parallelhelpers.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/matrixconverter.hh>
#include <dumux/common/typetraits/vector.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

#if DUMUX_HAVE_TRILINOS

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>
#include <Amesos2_MatrixAdapter_decl.hpp>

#define HAVE_MPI 1
#if HAVE_MPI
#include <mpi.h>
#endif

namespace Dumux {

namespace Detail {

// Provides GridView, DofMapper, dofCodim from LSTraits when canCommunicate=true,
// or inert dummy types when canCommunicate=false (SeqLinearSolverTraits).
template<class LSTraits, bool = LSTraits::canCommunicate>
struct AmesosGridTypes {
    struct GridView {};
    struct DofMapper { std::size_t size() const { return 0; } };
    static constexpr int dofCodim = 0;
};

template<class LSTraits>
struct AmesosGridTypes<LSTraits, true> {
    using GridView = typename LSTraits::GridView;
    using DofMapper = typename LSTraits::DofMapper;
    static constexpr int dofCodim = LSTraits::dofCodim;
};

// Helper to extract sub-vector types from MultiTypeBlockVector without instantiating
// std::tuple_element on non-tuple types.
template<class V, bool IsMT>
struct MTSubVecs { using Mom = V; using Mass = V; };
template<class V>
struct MTSubVecs<V, true> {
    using Mom = std::tuple_element_t<0, V>;
    using Mass = std::tuple_element_t<1, V>;
};

} // end namespace Detail

/*!
 * \ingroup Linear
 * \brief Direct linear solvers from Trilinos Amesos2
 *
 * Wraps MUMPS (or other Amesos2 backends) via Tpetra for parallel direct solves.
 * Handles both sequential and parallel (non-overlapping and overlapping) decompositions.
 * Also supports MultiTypeBlockMatrix/Vector (multidomain) via automatic conversion to scalar types.
 */
template<class LSTraits, class LATraits>
class DirectSolverAmesos2 : public LinearSolver
{
    using Matrix = typename LATraits::Matrix;
    using XVector = typename LATraits::Vector;
    using BVector = typename LATraits::Vector;

    using Scalar = double;
    using NodeType = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;
    using LocalIDType = Tpetra::Map<>::local_ordinal_type;
    using GlobalIDType = Tpetra::Map<>::global_ordinal_type;
    using TpetraMap = Tpetra::Map<LocalIDType, GlobalIDType, NodeType>;
    using SolverMatrix = Tpetra::CrsMatrix<Scalar, LocalIDType, GlobalIDType, NodeType>;
    using SolverVector = Tpetra::MultiVector<Scalar, LocalIDType, GlobalIDType, NodeType>;
    using AmeSolver = Amesos2::Solver<SolverMatrix, SolverVector>;
    using Graph = Tpetra::CrsGraph<LocalIDType, GlobalIDType, NodeType>;

    using GridTypes_ = Detail::AmesosGridTypes<LSTraits>;
    using GridView = typename GridTypes_::GridView;
    using DofMapper = typename GridTypes_::DofMapper;
    static constexpr int dofCodim = GridTypes_::dofCodim;

    // True when XVector is a MultiTypeBlockVector (multidomain assembler)
    static constexpr bool isMultiType = isMultiTypeBlockVector<XVector>::value;
    // Scalar DOFs per block DOF. MultiType is always converted to scalar (blockSize=1).
    static constexpr std::size_t blockSize = []() constexpr -> std::size_t {
        if constexpr (isMultiTypeBlockVector<XVector>::value) return 1;
        else return XVector::block_type::size();
    }();

    // Sub-vector types and scalar block sizes per subdomain (MultiType only).
    // When !isMultiType, Mom/Mass fall back to XVector — values are unused.
    using MomBlockVec_  = typename Detail::MTSubVecs<XVector, isMultiType>::Mom;
    using MassBlockVec_ = typename Detail::MTSubVecs<XVector, isMultiType>::Mass;
    static constexpr std::size_t momScalarBS_  = MomBlockVec_::block_type::size();
    static constexpr std::size_t massScalarBS_ = MassBlockVec_::block_type::size();

public:
    struct SolverResult : public Dune::InverseOperatorResult
    {
        SolverResult() = default;
        SolverResult(const SolverResult&) = default;
        SolverResult(SolverResult&&) = default;
        SolverResult(const Dune::InverseOperatorResult& o) : InverseOperatorResult(o) {}
        SolverResult(Dune::InverseOperatorResult&& o) : InverseOperatorResult(std::move(o)) {}
        operator bool() const { return this->converged; }
    };

    // Constructor for sequential (SeqLinearSolverTraits) or single-type solvers with no grid.
    // DOF maps are built lazily on the first solve() call.
    explicit DirectSolverAmesos2(const std::string& paramGroup = "")
    requires (!LSTraits::canCommunicate)
    : LinearSolver(paramGroup)
    , isParallel_(false)
    {
        comm_ = Teuchos::rcp(new Teuchos::MpiComm<int>(MPI_COMM_SELF));
        fos_ = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        rank_ = comm_->getRank();
        if (rank_ == 0)
            *fos_ << "\n" << Amesos2::version() << std::endl;
    }

    // Constructor for parallel/sequential solvers with a grid (canCommunicate=true).
    template<class GridGeometry>
    DirectSolverAmesos2(const GridGeometry& gridGeometry,
                        const GridView& gridView,
                        const DofMapper& dofMapper,
                        const std::string& paramGroup = "")
    requires (LSTraits::canCommunicate)
    : LinearSolver(paramGroup)
    , gridView_(gridView)
    , mapper_(&dofMapper)
    , isParallel_(gridView.comm().size() > 1)
    {
        comm_ = Teuchos::rcp(new Teuchos::MpiComm<int>(static_cast<MPI_Comm>(gridView.comm())));
        fos_ = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        rank_ = comm_->getRank();
        if (rank_ == 0)
            *fos_ << "\n" << Amesos2::version() << std::endl;

        buildGlobalDofIndices_(gridView, dofMapper);
        rowMap_ = buildRowMap_();
        colMap_ = buildColMap_();
    }

    // Constructor for parallel MultiType (multidomain) problems.
    // Takes both subdomain grid geometries to build combined parallel DOF maps.
    template<class MomGG, class MassGG>
    DirectSolverAmesos2(const MomGG& momGG, const MassGG& massGG,
                        const std::string& paramGroup = "")
    requires (isMultiType)
    : LinearSolver(paramGroup)
    , isParallel_(momGG.gridView().comm().size() > 1)
    {
        comm_ = Teuchos::rcp(new Teuchos::MpiComm<int>(
            static_cast<MPI_Comm>(momGG.gridView().comm())));
        fos_ = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        rank_ = comm_->getRank();
        if (rank_ == 0)
            *fos_ << "\n" << Amesos2::version() << std::endl;

        using MomTraits = LinearSolverTraits<MomGG>;
        using MassTraits = LinearSolverTraits<MassGG>;

        std::vector<GlobalIDType> momGIDs, massGIDs;
        std::vector<bool> momOwned, momGhost, massOwned, massGhost;

        const auto nMomGlobal = buildSubdomainGIDs_<MomTraits>(
            momGG.gridView(), momGG.dofMapper(), 0, momGIDs, momOwned, momGhost);
        const auto nMassGlobal = buildSubdomainGIDs_<MassTraits>(
            massGG.gridView(), massGG.dofMapper(), nMomGlobal, massGIDs, massOwned, massGhost);

        // Expand block-level GIDs to scalar level.
        // multiTypeToBCRSMatrix expands momentum DOF i to momScalarBS_ scalar rows (2 for 2D),
        // so localToGlobal_, isOwned_, isGhost_ must be sized at the scalar row count.
        nMomBlockLocal_  = momGIDs.size();
        nMassBlockLocal_ = massGIDs.size();
        nMomLocal_  = momScalarBS_  * nMomBlockLocal_;
        nMassLocal_ = massScalarBS_ * nMassBlockLocal_;

        const GlobalIDType nMomScalarGlobal  = static_cast<GlobalIDType>(momScalarBS_)  * nMomGlobal;
        const GlobalIDType nMassScalarGlobal = static_cast<GlobalIDType>(massScalarBS_) * nMassGlobal;
        numBlocksGlobal_ = static_cast<std::size_t>(nMomScalarGlobal + nMassScalarGlobal);

        localToGlobal_.resize(nMomLocal_ + nMassLocal_);
        isOwned_.resize(nMomLocal_ + nMassLocal_);
        isGhost_.resize(nMomLocal_ + nMassLocal_);

        for (std::size_t i = 0; i < nMomBlockLocal_; ++i)
            for (std::size_t k = 0; k < momScalarBS_; ++k) {
                localToGlobal_[momScalarBS_*i + k] = momScalarBS_ * momGIDs[i] + static_cast<GlobalIDType>(k);
                isOwned_[momScalarBS_*i + k] = momOwned[i];
                isGhost_[momScalarBS_*i + k] = momGhost[i];
            }

        for (std::size_t j = 0; j < nMassBlockLocal_; ++j)
            for (std::size_t k = 0; k < massScalarBS_; ++k) {
                // massGIDs[j] is the block GID offset from nMomGlobal; expand to scalar
                localToGlobal_[nMomLocal_ + massScalarBS_*j + k] =
                    nMomScalarGlobal + massScalarBS_ * (massGIDs[j] - nMomGlobal) + static_cast<GlobalIDType>(k);
                isOwned_[nMomLocal_ + massScalarBS_*j + k] = massOwned[j];
                isGhost_[nMomLocal_ + massScalarBS_*j + k] = massGhost[j];
            }

        rowMap_ = buildRowMap_();
        colMap_ = buildColMap_();

        if (isParallel_)
            setupParallelMultiTypeFuncs_<MomTraits, MassTraits>(momGG, massGG);
    }

    SolverResult solve(const Matrix& A, XVector& x, const BVector& b)
    { return solve_(A, x, b); }

    Scalar norm(const XVector& x) const
    {
        if constexpr (isMultiType) {
            if (!isParallel_ || localToGlobal_.empty())
            {
                auto y = VectorConverter<XVector>::multiTypeToBlockVector(x);
                return y.two_norm();
            }
            // Parallel: apply makeNonOverlappingConsistent to sum border DOF contributions
            // from all processes before computing the owned-DOF global norm.
            // Without this, border DOF partial contributions inflate the partial norm
            // and it stays constant even as the solution converges.
            MomBlockVec_ bMomBlock(nMomBlockLocal_);
            for (std::size_t i = 0; i < nMomBlockLocal_; ++i)
                bMomBlock[i] = x[Dune::Indices::_0][i];
            if (momRhsPreprocess_)
                momRhsPreprocess_(bMomBlock);

            MassBlockVec_ bMassBlock(nMassBlockLocal_);
            for (std::size_t j = 0; j < nMassBlockLocal_; ++j)
                bMassBlock[j] = x[Dune::Indices::_1][j];
            if (massRhsPreprocess_)
                massRhsPreprocess_(bMassBlock);

            double localSumSq = 0.0;
            for (std::size_t i = 0; i < nMomBlockLocal_; ++i)
                if (isOwned_[momScalarBS_*i])
                    for (std::size_t k = 0; k < momScalarBS_; ++k)
                    {
                        const auto v = bMomBlock[i][k];
                        localSumSq += v * v;
                    }

            // Mass owned contributions (preprocessed if needed)
            for (std::size_t j = 0; j < nMassBlockLocal_; ++j)
                if (isOwned_[nMomLocal_ + massScalarBS_*j])
                    for (std::size_t k = 0; k < massScalarBS_; ++k)
                    {
                        const auto v = bMassBlock[j][k];
                        localSumSq += v * v;
                    }

            double globalSumSq = 0.0;
            Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &localSumSq, &globalSumSq);
            using std::sqrt;
            return sqrt(globalSumSq);
        } else {
            if (!isParallel_ || localToGlobal_.empty())
                return x.two_norm();

            // Parallel non-MultiType: compute the true global residual norm.
            // The Newton solver expects the linear solver to perform all communication.
            // For non-overlapping decompositions: the assembled vector x contains only
            // each rank's partial contributions to shared border DOFs. We must first
            // accumulate those partial contributions from all ranks (making the vector
            // globally consistent), then compute the owned-DOF global norm.
            // For overlapping decompositions: owned DOFs already have complete values,
            // so only the MPI reduction is needed (no intra-norm communication).
            if constexpr (LSTraits::canCommunicate)
            {
                if (LSTraits::isNonOverlapping(gridView_))
                {
                    // Make a mutable copy and accumulate partial border-DOF contributions
                    // from all non-owning ranks onto the owners.
                    auto xcopy = x;
                    ParallelVectorHelper<GridView, DofMapper, dofCodim> vecHelper(gridView_, *mapper_);
                    if constexpr (requires { LSTraits::dofCodims; })
                        vecHelper.makeNonOverlappingConsistent(xcopy, LSTraits::dofCodims);
                    else
                        vecHelper.makeNonOverlappingConsistent(xcopy);

                    double localSumSq = 0.0;
                    for (std::size_t i = 0; i < xcopy.size(); ++i)
                        if (isOwned_[i])
                            for (std::size_t k = 0; k < blockSize; ++k)
                            {
                                const auto v = xcopy[i][k];
                                localSumSq += v * v;
                            }
                    double globalSumSq = 0.0;
                    Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &localSumSq, &globalSumSq);
                    using std::sqrt;
                    return sqrt(globalSumSq);
                }
            }

            // Overlapping or other parallel: owned-DOF global norm without extra communication
            double localSumSq = 0.0;
            for (std::size_t i = 0; i < x.size(); ++i)
                if (isOwned_[i])
                    for (std::size_t k = 0; k < blockSize; ++k)
                    {
                        const auto v = x[i][k];
                        localSumSq += v * v;
                    }
            double globalSumSq = 0.0;
            Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &localSumSq, &globalSumSq);
            using std::sqrt;
            return sqrt(globalSumSq);
        }
    }

    std::string name() const { return "Trilinos Amesos2 Direct Solver"; }

private:
    /*!
     * Compute localToGlobal_[i] = Tpetra global block index for Dune local block DOF i.
     * Scalar Tpetra GID for (block i, sub-index k) = localToGlobal_[i] * blockSize + k.
     * Also populates isOwned_[i] and isGhost_[i].
     */
    void buildGlobalDofIndices_(const GridView& gridView, const DofMapper& mapper)
    {
        const std::size_t nLocal = mapper.size();
        localToGlobal_.resize(nLocal, std::numeric_limits<GlobalIDType>::max());
        isOwned_.assign(nLocal, false);
        isGhost_.assign(nLocal, false);

        if constexpr (LSTraits::canCommunicate)
        {
            if (isParallel_)
            {
                // Use ParallelISTLHelper to determine ownership and ghost status
                ParallelISTLHelper<LSTraits> pHelper(gridView, mapper);
                for (std::size_t i = 0; i < nLocal; ++i)
                {
                    isOwned_[i] = pHelper.isOwned(i);
                    isGhost_[i] = pHelper.isGhost(i);
                }

                // Count owned DOFs, gather counts, compute global offset
                GlobalIDType numOwned = static_cast<GlobalIDType>(
                    std::count(isOwned_.begin(), isOwned_.end(), true));
                const int commSize = gridView.comm().size();
                const int myRank = gridView.comm().rank();
                std::vector<GlobalIDType> counts(commSize);
                gridView.comm().allgather(&numOwned, 1, counts.data());

                GlobalIDType offset = 0;
                for (int p = 0; p < myRank; ++p)
                    offset += counts[p];
                numBlocksGlobal_ = 0;
                for (auto c : counts) numBlocksGlobal_ += c;

                // Assign consecutive global block IDs to owned DOFs
                for (std::size_t i = 0; i < nLocal; ++i)
                    if (isOwned_[i])
                        localToGlobal_[i] = offset++;

                // Propagate global IDs from owners to ghost DOFs via min communication
                // (ghosts start with max value, owner sends real ID, min picks it up)
                if constexpr (requires { LSTraits::dofCodims; }) {
                        MultiCodimVectorCommDataHandleMin<DofMapper, std::vector<GlobalIDType>, GridView::dimension, GlobalIDType>
                        handle(mapper, localToGlobal_, LSTraits::dofCodims);
                    // InteriorBorder_All_Interface: owners (interior/border) send IDs to ALL
                    // (including ghost). Using All_All_Interface causes owned DOF global IDs
                    // to propagate back to ghost DOFs that are ALUGrid's duplicate representation
                    // of the same physical border entity (appearing as both BorderEntity and
                    // GhostEntity on the same rank due to ghost element sub-entity indexing).
                    gridView.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
                } else {
                    VectorCommDataHandleMin<DofMapper, std::vector<GlobalIDType>, dofCodim, GlobalIDType>
                        handle(mapper, localToGlobal_);
                    // InteriorBorder_All_Interface: owners (interior/border) send IDs to ALL
                    // (including ghost). Using All_All_Interface causes owned DOF global IDs
                    // to propagate back to ghost DOFs that are ALUGrid's duplicate representation
                    // of the same physical border entity (appearing as both BorderEntity and
                    // GhostEntity on the same rank due to ghost element sub-entity indexing).
                    gridView.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
                }

                return;
            }
        }

        // Sequential: all DOFs are owned with consecutive IDs
        numBlocksGlobal_ = nLocal;
        std::fill(isOwned_.begin(), isOwned_.end(), true);
        for (std::size_t i = 0; i < nLocal; ++i)
            localToGlobal_[i] = static_cast<GlobalIDType>(i);
    }

    // Initialize DOF maps for sequential use from a known scalar DOF count.
    // Called lazily on first solve for !canCommunicate solvers.
    void initSeqMaps_(std::size_t nDofs)
    {
        numBlocksGlobal_ = nDofs;
        localToGlobal_.resize(nDofs);
        isOwned_.assign(nDofs, true);
        isGhost_.assign(nDofs, false);
        for (std::size_t i = 0; i < nDofs; ++i)
            localToGlobal_[i] = static_cast<GlobalIDType>(i);
        rowMap_ = buildRowMap_();
        colMap_ = buildColMap_();
    }

    // Build global DOF indices for one subdomain.
    // Returns the total number of global block DOFs in this subdomain.
    // gids[i] = global block ID for local DOF i; offset makes IDs globally unique.
    template<class SubLSTraits, class GV, class Mapper>
    GlobalIDType buildSubdomainGIDs_(const GV& gridView, const Mapper& mapper,
                                     GlobalIDType globalOffset,
                                     std::vector<GlobalIDType>& gids,
                                     std::vector<bool>& owned,
                                     std::vector<bool>& ghost)
    {
        const std::size_t nLocal = mapper.size();
        gids.resize(nLocal, std::numeric_limits<GlobalIDType>::max());
        owned.assign(nLocal, false);
        ghost.assign(nLocal, false);

        if constexpr (SubLSTraits::canCommunicate)
        {
            if (isParallel_)
            {
                ParallelISTLHelper<SubLSTraits> pHelper(gridView, mapper);
                for (std::size_t i = 0; i < nLocal; ++i)
                {
                    owned[i] = pHelper.isOwned(i);
                    ghost[i] = pHelper.isGhost(i);
                }

                GlobalIDType numOwned = static_cast<GlobalIDType>(
                    std::count(owned.begin(), owned.end(), true));
                const int commSize = gridView.comm().size();
                const int myRank = gridView.comm().rank();
                std::vector<GlobalIDType> counts(commSize);
                gridView.comm().allgather(&numOwned, 1, counts.data());

                GlobalIDType localOffset = globalOffset;
                for (int p = 0; p < myRank; ++p)
                    localOffset += counts[p];

                GlobalIDType nGlobal = 0;
                for (auto c : counts) nGlobal += c;

                for (std::size_t i = 0; i < nLocal; ++i)
                    if (owned[i])
                        gids[i] = localOffset++;

                // Propagate owned GIDs to ghost DOFs via min communication
                if constexpr (requires { SubLSTraits::dofCodims; }) {
                    MultiCodimVectorCommDataHandleMin<Mapper, std::vector<GlobalIDType>, GV::dimension, GlobalIDType>
                        handle(mapper, gids, SubLSTraits::dofCodims);
                    // InteriorBorder_All_Interface: owners (interior/border) send IDs to ALL
                    // (including ghost). Using All_All_Interface causes owned DOF global IDs
                    // to propagate back to ghost DOFs that are ALUGrid's duplicate representation
                    // of the same physical border entity (appearing as both BorderEntity and
                    // GhostEntity on the same rank due to ghost element sub-entity indexing).
                    gridView.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
                } else {
                    constexpr int codim = SubLSTraits::dofCodim;
                    VectorCommDataHandleMin<Mapper, std::vector<GlobalIDType>, codim, GlobalIDType>
                        handle(mapper, gids);
                    // InteriorBorder_All_Interface: owners (interior/border) send IDs to ALL
                    // (including ghost). Using All_All_Interface causes owned DOF global IDs
                    // to propagate back to ghost DOFs that are ALUGrid's duplicate representation
                    // of the same physical border entity (appearing as both BorderEntity and
                    // GhostEntity on the same rank due to ghost element sub-entity indexing).
                    gridView.communicate(handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
                }
                return nGlobal;
            }
        }

        // Sequential path: all DOFs owned with consecutive IDs
        std::fill(owned.begin(), owned.end(), true);
        for (std::size_t i = 0; i < nLocal; ++i)
            gids[i] = globalOffset + static_cast<GlobalIDType>(i);
        return static_cast<GlobalIDType>(nLocal);
    }

    // Set up block-level RHS preprocessing and ghost broadcast closures for both subdomains.
    // Matrix border-row accumulation is handled by Tpetra global assembly (convertMatrixGlobal_).
    template<class MomTraits, class MassTraits, class MomGG, class MassGG>
    void setupParallelMultiTypeFuncs_(const MomGG& momGG, const MassGG& massGG)
    {
        using MomGVType = typename MomGG::GridView;
        using MomMapperType = typename MomGG::DofMapper;
        using MassMapperType = typename MassGG::DofMapper;

        const auto momGV = momGG.gridView();
        const auto* momMapper = &momGG.dofMapper();
        const auto massGV = massGG.gridView();
        const auto* massMapper = &massGG.dofMapper();

        // Block-level momentum ghost broadcast (operates on MomBlockVec_ = BlockVector<FV<dim>>)
        if constexpr (requires { MomTraits::dofCodims; }) {
            constexpr int numMomCodims = MomGVType::dimension;
            momGhostBroadcastBlock_ = [momGV, momMapper](MomBlockVec_& vec) {
                MultiCodimVectorCommDataHandleSum<MomMapperType, MomBlockVec_, numMomCodims>
                    handle(*momMapper, vec, MomTraits::dofCodims);
                momGV.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
            };
        } else {
            constexpr int momCodim = MomTraits::dofCodim;
            momGhostBroadcastBlock_ = [momGV, momMapper](MomBlockVec_& vec) {
                VectorCommDataHandleSum<MomMapperType, MomBlockVec_, momCodim>
                    handle(*momMapper, vec);
                momGV.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
            };
        }

        // Block-level mass ghost broadcast
        constexpr int massCodim = MassTraits::dofCodim;
        using MassGVType = typename MassGG::GridView;
        massGhostBroadcastBlock_ = [massGV, massMapper](MassBlockVec_& vec) {
            VectorCommDataHandleSum<MassMapperType, MassBlockVec_, massCodim>
                handle(*massMapper, vec);
            massGV.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
        };

        // For vertex-based mass DOFs (Box, codim=dim) on non-overlapping decompositions,
        // CVFELocalAssembler skips ghost-element contributions to border vertex DOFs,
        // so owned border mass rows are INCOMPLETE — same as momentum.
        // For face-based mass (FCDiamond, codim=1), interior assembly captures both sides.
        constexpr bool massAtVertices = (MassTraits::dofCodim == MassGVType::dimension);
        massRowsAreComplete_ = !(massAtVertices && MassTraits::isNonOverlapping(massGV));

        if (!massRowsAreComplete_) {
            // RHS: sum border mass DOF contributions from all processes
            massRhsPreprocess_ = [massGV, massMapper](MassBlockVec_& bMass) {
                ParallelVectorHelper<MassGVType, MassMapperType, massCodim>
                    vecHelper(massGV, *massMapper);
                vecHelper.makeNonOverlappingConsistent(bMass);
            };
            // Extra-count communication for pre-allocating Tpetra capacity
            massExtraCountComm_ = [massGV, massMapper](std::vector<std::size_t>& counts) {
                VectorCommDataHandleSum<MassMapperType, std::vector<std::size_t>, massCodim, std::size_t>
                    handle(*massMapper, counts);
                massGV.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
            };
        }

        // Block-level momentum RHS preprocessing for non-overlapping decomposition:
        // sum border vertex DOF contributions across processes before scalar conversion.
        if (MomTraits::isNonOverlapping(momGV)) {
            if constexpr (requires { MomTraits::dofCodims; }) {
                momRhsPreprocess_ = [momGV, momMapper](MomBlockVec_& bMom) {
                    ParallelVectorHelper<MomGVType, MomMapperType, MomTraits::dofCodim>
                        vecHelper(momGV, *momMapper);
                    vecHelper.makeNonOverlappingConsistent(bMom, MomTraits::dofCodims);
                };
            } else {
                constexpr int momCodim = MomTraits::dofCodim;
                momRhsPreprocess_ = [momGV, momMapper](MomBlockVec_& bMom) {
                    ParallelVectorHelper<MomGVType, MomMapperType, momCodim>
                        vecHelper(momGV, *momMapper);
                    vecHelper.makeNonOverlappingConsistent(bMom);
                };
            }
        }

        // Communication closure to sum up extra entry counts for border momentum DOFs.
        // Non-owning processes send their partial row sizes; the owner accumulates the total
        // so it can pre-allocate enough capacity for global Tpetra assembly.
        if constexpr (requires { MomTraits::dofCodims; }) {
            constexpr int numMomCodims = MomGVType::dimension;
            momExtraCountComm_ = [momGV, momMapper](std::vector<std::size_t>& counts) {
                MultiCodimVectorCommDataHandleSum<MomMapperType, std::vector<std::size_t>, numMomCodims, std::size_t>
                    handle(*momMapper, counts, MomTraits::dofCodims);
                momGV.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
            };
        } else {
            constexpr int momCodim = MomTraits::dofCodim;
            momExtraCountComm_ = [momGV, momMapper](std::vector<std::size_t>& counts) {
                VectorCommDataHandleSum<MomMapperType, std::vector<std::size_t>, momCodim, std::size_t>
                    handle(*momMapper, counts);
                momGV.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
            };
        }
    }

    /*!
     * Row map: owned scalar DOFs only.
     * Scalar GID for (block i, sub k) = localToGlobal_[i] * blockSize + k.
     */
    Teuchos::RCP<const TpetraMap> buildRowMap_()
    {
        const std::size_t nLocal = localToGlobal_.size();
        std::vector<GlobalIDType> ownedGIDs;
        ownedGIDs.reserve(
            static_cast<std::size_t>(std::count(isOwned_.begin(), isOwned_.end(), true)) * blockSize);
        for (std::size_t i = 0; i < nLocal; ++i)
            if (isOwned_[i])
                for (std::size_t k = 0; k < blockSize; ++k)
                    ownedGIDs.push_back(localToGlobal_[i] * blockSize + k);

        const auto numGlobalScalarRows =
            static_cast<Tpetra::global_size_t>(numBlocksGlobal_ * blockSize);
        Teuchos::ArrayView<const GlobalIDType> gidsView(ownedGIDs.data(), ownedGIDs.size());
        return Teuchos::rcp(new TpetraMap(numGlobalScalarRows, gidsView, GlobalIDType{0}, comm_));
    }

    /*!
     * Col map: ALL local Dune DOFs (owned + ghost), Dune-local-index order.
     * Invariant: colMap_.getLocalElement(localToGlobal_[i]*blockSize + k) == i*blockSize + k
     * This lets us use Dune local block index directly as Tpetra local column block index.
     */
    Teuchos::RCP<const TpetraMap> buildColMap_()
    {
        const std::size_t nLocal = localToGlobal_.size();
        std::vector<GlobalIDType> allGIDs;
        allGIDs.reserve(nLocal * blockSize);
        for (std::size_t i = 0; i < nLocal; ++i)
            for (std::size_t k = 0; k < blockSize; ++k)
                allGIDs.push_back(localToGlobal_[i] * blockSize + k);

        Teuchos::ArrayView<const GlobalIDType> gidsView(allGIDs.data(), allGIDs.size());
        return Teuchos::rcp(new TpetraMap(
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
            gidsView, GlobalIDType{0}, comm_));
    }

    /*!
     * Build Tpetra CrsGraph from the actual Dune matrix (which may have been
     * extended by extendMatrix for parallel runs).
     * Templated on AMatrix to handle both plain BCRSMatrix and scalar-converted types.
     */
    /*!
     * Global assembly for non-overlapping parallel direct solves.
     *
     * Submits every non-ghost block row (owned rows + non-owned non-ghost border rows)
     * using Tpetra insertGlobalValues.  Non-owned rows become non-local entries that
     * fillComplete redistributes to the owning rank via global assembly (ADD).
     *
     * Correctness: DuMux ghost elements do NOT contribute physical stiffness to non-ghost
     * DOF rows (cvfelocalassembler.hh skips `assembleJacobianAndResidualImpl` for ghost
     * elements).  Therefore each element's stiffness is provided by exactly one rank
     * (its owner), with no double-counting:
     *   Rank A (owner of border DOFs): K_AA + K_AB   (own real elements only)
     *   Rank B (non-local inserts):    K_BA + K_BB_v  (all rank-B real elements adjacent
     *                                                    to the border DOFs, including
     *                                                    vertex-only adjacent ones missing
     *                                                    from the 1-layer ghost)
     *   After fillComplete ADD:  K_AA+K_AB+K_BA+K_BB_v = K_total  ✓
     */
    template<class AMatrix>
    Teuchos::RCP<SolverMatrix> convertMatrixGlobalNonOverlapping_(const AMatrix& A)
    {
        using BlockType = typename AMatrix::block_type;
        static constexpr int BS_I = BlockType::rows;
        static constexpr int BS_J = BlockType::cols;

        // Pass 1: count owned-row entries and communicate non-owned row sizes to owners.
        const std::size_t nTpetraRows = rowMap_->getLocalNumElements();
        Teuchos::ArrayRCP<std::size_t> numEntPerRow(nTpetraRows, std::size_t{0});

        // Extra entry counts that non-owned border rows will contribute to owners.
        std::vector<std::size_t> extraCounts(localToGlobal_.size(), 0);

        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            const std::size_t r = rowIt.index();
            if (isGhost_[r]) continue;           // skip ghost rows entirely
            const std::size_t nBlockCols = rowIt->size();

            if (isOwned_[r])
            {
                for (int ki = 0; ki < BS_I; ++ki)
                {
                    const LocalIDType tRow = rowMap_->getLocalElement(
                        static_cast<GlobalIDType>(localToGlobal_[r] * BS_I + ki));
                    numEntPerRow[tRow] = nBlockCols * BS_J;
                }
            }
            else
            {
                // Non-owned non-ghost border row → will become a non-local insert.
                // Tell the owner how many entries to expect.
                extraCounts[r] = nBlockCols; // block-column count
            }
        }

        // Communicate extra counts: non-owners → owners (border entities).
        if constexpr (requires { LSTraits::dofCodims; })
        {
            constexpr auto dofCodims = LSTraits::dofCodims;
            MultiCodimVectorCommDataHandleSum<DofMapper, std::vector<std::size_t>,
                                              GridView::dimension, std::size_t>
                handle(*mapper_, extraCounts, dofCodims);
            gridView_.communicate(handle, Dune::InteriorBorder_InteriorBorder_Interface,
                                  Dune::ForwardCommunication);
        }
        else
        {
            VectorCommDataHandleSum<DofMapper, std::vector<std::size_t>, dofCodim, std::size_t>
                handle(*mapper_, extraCounts);
            gridView_.communicate(handle, Dune::InteriorBorder_InteriorBorder_Interface,
                                  Dune::ForwardCommunication);
        }

        for (std::size_t r = 0; r < localToGlobal_.size(); ++r)
            if (extraCounts[r] > 0 && isOwned_[r])
                for (int ki = 0; ki < BS_I; ++ki)
                {
                    const LocalIDType tRow = rowMap_->getLocalElement(
                        static_cast<GlobalIDType>(localToGlobal_[r] * BS_I + ki));
                    numEntPerRow[tRow] += extraCounts[r] * BS_J;
                }

        // Pass 2: insert values using global row and column indices.
        auto AA = Teuchos::rcp(new SolverMatrix(rowMap_, numEntPerRow()));

        Teuchos::Array<GlobalIDType> gCols;
        Teuchos::Array<Scalar> vals;

        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            const std::size_t r = rowIt.index();
            if (isGhost_[r]) continue;

            const GlobalIDType gBlockRow = static_cast<GlobalIDType>(localToGlobal_[r]);

            for (int ki = 0; ki < BS_I; ++ki)
            {
                const GlobalIDType globalRow = gBlockRow * BS_I + ki;
                gCols.clear();
                vals.clear();
                for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                {
                    const GlobalIDType gBlockCol =
                        static_cast<GlobalIDType>(localToGlobal_[colIt.index()]);
                    for (int kj = 0; kj < BS_J; ++kj)
                    {
                        gCols.push_back(gBlockCol * BS_J + kj);
                        vals.push_back((*colIt)[ki][kj]);
                    }
                }
                AA->insertGlobalValues(globalRow, gCols(), vals());
            }
        }

        AA->fillComplete(rowMap_, rowMap_);
        return AA;
    }

    template<class AMatrix>
    Teuchos::RCP<const Graph> buildGraph_(const AMatrix& A)
    {
        const std::size_t nLocalTpetraRows = rowMap_->getLocalNumElements();

        // Count scalar entries per Tpetra row
        Teuchos::ArrayRCP<std::size_t> numEntPerRow(nLocalTpetraRows, 0);
        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            const std::size_t blockRow = rowIt.index();
            if (!isOwned_[blockRow]) continue;
            const std::size_t nBlockCols = rowIt->size();
            for (std::size_t ki = 0; ki < blockSize; ++ki)
            {
                const LocalIDType tRow = rowMap_->getLocalElement(
                    static_cast<GlobalIDType>(localToGlobal_[blockRow] * blockSize + ki));
                numEntPerRow[tRow] = nBlockCols * blockSize;
            }
        }

        // Build graph using GLOBAL column indices so that ghost-column entries
        // (border DOF rows with ghost interior DOF columns from ghost element assembly)
        // are included. Use dynamic profile (no pre-allocation) to avoid any static
        // capacity limits silently dropping ghost-column entries.
        auto graph = Teuchos::rcp(new Graph(rowMap_, numEntPerRow()));

        Teuchos::Array<GlobalIDType> globalColIndices;
        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            const std::size_t blockRow = rowIt.index();
            if (!isOwned_[blockRow]) continue;

            globalColIndices.clear();
            globalColIndices.reserve(rowIt->size() * blockSize);
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                for (std::size_t kj = 0; kj < blockSize; ++kj)
                    globalColIndices.push_back(
                        static_cast<GlobalIDType>(localToGlobal_[colIt.index()] * blockSize + kj));

            for (std::size_t ki = 0; ki < blockSize; ++ki)
            {
                const GlobalIDType gRow =
                    static_cast<GlobalIDType>(localToGlobal_[blockRow] * blockSize + ki);
                graph->insertGlobalIndices(gRow, globalColIndices());
            }
        }

        graph->fillComplete(rowMap_, rowMap_);
        return graph;
    }

    SolverResult solve_(const Matrix& A, XVector& x, const BVector& b)
    {
        if constexpr (isMultiType)
        {
            // Flatten MultiTypeBlockMatrix/Vector to scalar BCRSMatrix/BlockVector for Tpetra
            auto Ascalar = MatrixConverter<Matrix>::multiTypeToBCRSMatrix(A);
            auto bscalar = VectorConverter<BVector>::multiTypeToBlockVector(b);

            // Build DOF maps lazily on first call for the sequential (no-arg constructor) path
            if (localToGlobal_.empty())
                initSeqMaps_(Ascalar.N());

            // Block-level RHS preprocessing: sum border DOF contributions across processes.
            if (isParallel_ && momRhsPreprocess_)
            {
                MomBlockVec_ bMomBlock(nMomBlockLocal_);
                for (std::size_t i = 0; i < nMomBlockLocal_; ++i)
                    bMomBlock[i] = b[Dune::Indices::_0][i];
                momRhsPreprocess_(bMomBlock);
                for (std::size_t i = 0; i < nMomBlockLocal_; ++i)
                    for (std::size_t k = 0; k < momScalarBS_; ++k)
                        bscalar[momScalarBS_*i + k] = bMomBlock[i][k];
            }
            // For vertex-based mass (Box): also sum border mass DOF contributions
            if (isParallel_ && massRhsPreprocess_)
            {
                MassBlockVec_ bMassBlock(nMassBlockLocal_);
                for (std::size_t j = 0; j < nMassBlockLocal_; ++j)
                    bMassBlock[j] = b[Dune::Indices::_1][j];
                massRhsPreprocess_(bMassBlock);
                for (std::size_t j = 0; j < nMassBlockLocal_; ++j)
                    for (std::size_t k = 0; k < massScalarBS_; ++k)
                        bscalar[nMomLocal_ + massScalarBS_*j + k] = bMassBlock[j][k];
            }

            // Build Tpetra matrix. Parallel: use global assembly so that non-local border-vertex
            // row contributions are accumulated by fillComplete (replaces extendMatrix+sumEntries).
            Teuchos::RCP<SolverMatrix> AA;
            if (isParallel_)
                AA = convertMatrixGlobal_(Ascalar);
            else
            {
                auto graph = buildGraph_(Ascalar);
                AA = convertMatrix_(Ascalar, graph);
            }

            auto BB = convertBVector_(bscalar);
            auto XX = convertXVector_();

            Teuchos::RCP<AmeSolver> solver = Amesos2::create<SolverMatrix, SolverVector>("MUMPS", AA, XX, BB);
            solver->setParameters(Teuchos::rcpFromRef(params_));
            solver->symbolicFactorization().numericFactorization().solve();

            auto xscalar = bscalar;
            retrieveXVector_(XX, xscalar);

            // Ghost broadcast at block level: split scalar solution into block sub-vectors,
            // zero non-owned DOFs, communicate, then reassemble to scalar for VectorConverter.
            if (isParallel_)
            {
                MomBlockVec_  xMomBlock(nMomBlockLocal_);
                MassBlockVec_ xMassBlock(nMassBlockLocal_);

                for (std::size_t i = 0; i < nMomBlockLocal_; ++i)
                    for (std::size_t k = 0; k < momScalarBS_; ++k)
                        xMomBlock[i][k] = xscalar[momScalarBS_*i + k];
                for (std::size_t j = 0; j < nMassBlockLocal_; ++j)
                    for (std::size_t k = 0; k < massScalarBS_; ++k)
                        xMassBlock[j][k] = xscalar[nMomLocal_ + massScalarBS_*j + k];

                for (std::size_t i = 0; i < nMomBlockLocal_; ++i)
                    if (!isOwned_[momScalarBS_*i]) xMomBlock[i] = 0;
                for (std::size_t j = 0; j < nMassBlockLocal_; ++j)
                    if (!isOwned_[nMomLocal_ + massScalarBS_*j]) xMassBlock[j] = 0;

                momGhostBroadcastBlock_(xMomBlock);
                massGhostBroadcastBlock_(xMassBlock);

                for (std::size_t i = 0; i < nMomBlockLocal_; ++i)
                    for (std::size_t k = 0; k < momScalarBS_; ++k)
                        xscalar[momScalarBS_*i + k] = xMomBlock[i][k];
                for (std::size_t j = 0; j < nMassBlockLocal_; ++j)
                    for (std::size_t k = 0; k < massScalarBS_; ++k)
                        xscalar[nMomLocal_ + massScalarBS_*j + k] = xMassBlock[j][k];
            }

            VectorConverter<XVector>::retrieveValues(x, xscalar);
        }
        else
        {
            auto Acopy = A;
            auto bcopy = b;

            if constexpr (LSTraits::canCommunicate)
            {
                if (isParallel_)
                {
                    if (LSTraits::isNonOverlapping(gridView_))
                    {
                        // For the non-overlapping case with MUMPS/Amesos2 direct solver,
                        // use Tpetra global assembly: each rank submits its non-ghost rows
                        // (owned rows + non-owned non-ghost border rows as non-local inserts).
                        // fillComplete redistributes non-local rows to their owners.
                        //
                        // Why no double-counting: ghost element physical stiffness is NOT
                        // assembled by the DuMux assembler (ghost elements only set identity
                        // for ghost DOF rows, see cvfelocalassembler.hh). So:
                        //   Rank A submits: K_AA + K_AB (own real elements only, no ghost)
                        //   Rank B submits: K_BA + K_BB_vertex (non-local for border DOF rows)
                        //   After fillComplete ADD: K_AA+K_AB+K_BA+K_BB_vertex = K_total ✓
                        // This handles both face-adjacent AND vertex-only-adjacent elements
                        // that would be missed by a 1-layer face-based ghost layer.
                        //
                        // RHS: the non-owned border rows have partial residuals from rank B
                        // that must also be accumulated onto the owner.
                        ParallelVectorHelper<GridView, DofMapper, dofCodim> vecHelper(gridView_, *mapper_);
                        if constexpr (requires { LSTraits::dofCodims; }) {
                            constexpr auto dofCodims = LSTraits::dofCodims;
                            vecHelper.makeNonOverlappingConsistent(bcopy, dofCodims);
                        } else {
                            vecHelper.makeNonOverlappingConsistent(bcopy);
                        }
                    }
                    // For overlapping decompositions (e.g. BOX with YaspGrid overlap=1):
                    // The overlap layer already provides each process with a complete local
                    // neighbourhood for all interior/border DOFs — no preprocessing needed.
                    // Non-owned (overlap/ghost) DOF rows are simply skipped during Tpetra conversion.
                }
            }

            // Use global assembly for parallel non-overlapping: non-ghost rows from rank B
            // are submitted as non-local inserts and redistributed by fillComplete.
            // For sequential or overlapping, use the graph-based path.
            Teuchos::RCP<SolverMatrix> AA;
            if (isParallel_ && LSTraits::canCommunicate && LSTraits::isNonOverlapping(gridView_))
                AA = convertMatrixGlobalNonOverlapping_(Acopy);
            else
            {
                auto graph = buildGraph_(Acopy);
                AA = convertMatrix_(Acopy, graph);
            }
            auto BB = convertBVector_(bcopy);
            auto XX = convertXVector_();

            Teuchos::RCP<AmeSolver> solver = Amesos2::create<SolverMatrix, SolverVector>("MUMPS", AA, XX, BB);
            solver->setParameters(Teuchos::rcpFromRef(params_));
            solver->symbolicFactorization().numericFactorization().solve();

            retrieveXVector_(XX, x);

            // Distribute solution: zero non-owned entries, then sum-communicate.
            // Non-owners send 0, owner sends the real value; each DOF ends up with the owner's value.
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
        }

        SolverResult result;
        result.converged = true;
        return result;
    }

    template<class AMatrix>
    Teuchos::RCP<SolverMatrix> convertMatrix_(const AMatrix& A, const Teuchos::RCP<const Graph>& graph)
    {
        Teuchos::RCP<SolverMatrix> AA = Teuchos::rcp(new SolverMatrix(graph));

        using BlockType = typename AMatrix::block_type;
        static constexpr int BS_I = BlockType::rows;
        static constexpr int BS_J = BlockType::cols;

        // Use global column indices (matching the graph which was built with insertGlobalIndices).
        // replaceGlobalValues works on a fillComplete'd matrix built from a global-index graph.
        Teuchos::Array<GlobalIDType> gCols;
        Teuchos::Array<Scalar> vals;

        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            const std::size_t blockRow = rowIt.index();
            if (!isOwned_[blockRow]) continue;

            const GlobalIDType gBlockRow = static_cast<GlobalIDType>(localToGlobal_[blockRow]);

            for (int ki = 0; ki < BS_I; ++ki)
            {
                const GlobalIDType globalRow = gBlockRow * BS_I + ki;
                gCols.clear();
                vals.clear();
                for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                {
                    const GlobalIDType gBlockCol = static_cast<GlobalIDType>(localToGlobal_[colIt.index()]);
                    for (int kj = 0; kj < BS_J; ++kj)
                    {
                        gCols.push_back(gBlockCol * BS_J + kj);
                        vals.push_back((*colIt)[ki][kj]);
                    }
                }
                // Use sumIntoGlobalValues: multiple local DOF columns may map to the
                // same global column GID when ALUGrid represents the same physical
                // border entity as both a BorderEntity and a GhostEntity sub-entity
                // of a ghost element. sumInto correctly accumulates both contributions.
                AA->sumIntoGlobalValues(globalRow, gCols(), vals());
            }
        }

        AA->fillComplete();

        return AA;
    }

    /*!
     * Build and fill Tpetra matrix using two-pass global assembly for parallel MultiType.
     * Pass 1: compute per-row Tpetra capacity, communicating non-local border row sizes.
     * Pass 2: insert entries using insertGlobalValues; fillComplete handles communication.
     * Momentum rows: submit if !isGhost_ (owned + border-non-owner both contribute).
     * Mass rows: submit if !isGhost_ when rows are incomplete (vertex-based, e.g. Box),
     *            or only if isOwned_ when rows are complete (face-based, e.g. FCDiamond).
     * AMatrix must be the scalar BCRSMatrix<FM<1,1>> from multiTypeToBCRSMatrix.
     */
    template<class AMatrix>
    Teuchos::RCP<SolverMatrix> convertMatrixGlobal_(const AMatrix& A)
    {
        // Pass 1: Per-row entry count + communication of non-local border row sizes.
        const std::size_t nTpetraRows = rowMap_->getLocalNumElements();
        Teuchos::ArrayRCP<std::size_t> numEntPerRow(nTpetraRows, std::size_t{0});

        std::vector<std::size_t> momExtraCounts(nMomBlockLocal_, 0);
        std::vector<std::size_t> massExtraCounts(massRowsAreComplete_ ? 0 : nMassBlockLocal_, 0);

        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            const std::size_t r = rowIt.index();
            const bool isMomRow = r < nMomLocal_;
            if (isMomRow)
            {
                if (isOwned_[r])
                {
                    numEntPerRow[rowMap_->getLocalElement(localToGlobal_[r])] = rowIt->size();
                }
                else if (!isGhost_[r] && r % momScalarBS_ == 0)
                {
                    momExtraCounts[r / momScalarBS_] = rowIt->size();
                }
            }
            else
            {
                const std::size_t massR = r - nMomLocal_;
                if (isOwned_[r])
                {
                    numEntPerRow[rowMap_->getLocalElement(localToGlobal_[r])] = rowIt->size();
                }
                else if (!massRowsAreComplete_ && !isGhost_[r] && massR % massScalarBS_ == 0)
                {
                    massExtraCounts[massR / massScalarBS_] = rowIt->size();
                }
            }
        }

        if (momExtraCountComm_)  momExtraCountComm_(momExtraCounts);
        if (massExtraCountComm_) massExtraCountComm_(massExtraCounts);

        for (std::size_t i = 0; i < nMomBlockLocal_; ++i)
            if (momExtraCounts[i] > 0)
                for (std::size_t k = 0; k < momScalarBS_; ++k)
                    if (isOwned_[momScalarBS_*i + k])
                        numEntPerRow[rowMap_->getLocalElement(localToGlobal_[momScalarBS_*i + k])]
                            += momExtraCounts[i];

        for (std::size_t j = 0; j < massExtraCounts.size(); ++j)
            if (massExtraCounts[j] > 0)
                for (std::size_t k = 0; k < massScalarBS_; ++k)
                    if (isOwned_[nMomLocal_ + massScalarBS_*j + k])
                        numEntPerRow[rowMap_->getLocalElement(
                            localToGlobal_[nMomLocal_ + massScalarBS_*j + k])]
                            += massExtraCounts[j];

        // Pass 2: Create matrix with correct capacity and insert values
        auto AA = Teuchos::rcp(new SolverMatrix(rowMap_, numEntPerRow()));

        Teuchos::Array<GlobalIDType> gCols;
        Teuchos::Array<Scalar> vals;

        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            const std::size_t r = rowIt.index();
            const bool isMomRow = r < nMomLocal_;
            bool shouldSubmit;
            if (isMomRow)
                shouldSubmit = !isGhost_[r];
            else
                shouldSubmit = massRowsAreComplete_ ? isOwned_[r] : !isGhost_[r];
            if (!shouldSubmit) continue;

            const GlobalIDType gRow = localToGlobal_[r];
            gCols.clear();
            vals.clear();
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
            {
                gCols.push_back(localToGlobal_[colIt.index()]);
                vals.push_back((*colIt)[0][0]);
            }
            AA->insertGlobalValues(gRow, gCols(), vals());
        }

        AA->fillComplete(rowMap_, rowMap_);
        return AA;
    }

    template<class BVec>
    Teuchos::RCP<SolverVector> convertBVector_(const BVec& b)
    {
        Teuchos::RCP<SolverVector> BB = Teuchos::rcp(new SolverVector(rowMap_, 1, true));
        for (std::size_t blockIdx = 0; blockIdx < b.size(); ++blockIdx)
        {
            if (!isOwned_[blockIdx]) continue;
            for (std::size_t j = 0; j < blockSize; ++j)
            {
                const LocalIDType tRow = rowMap_->getLocalElement(
                    static_cast<GlobalIDType>(localToGlobal_[blockIdx] * blockSize + j));
                BB->replaceLocalValue(tRow, 0, b[blockIdx][j]);
            }
        }
        return BB;
    }

    Teuchos::RCP<SolverVector> convertXVector_()
    {
        return Teuchos::rcp(new SolverVector(rowMap_, 1, true));
    }

    template<class XVec>
    void retrieveXVector_(const Teuchos::RCP<SolverVector>& XX, XVec& x)
    {
        const auto view = XX->getLocalViewHost(Tpetra::Access::ReadOnly);
        const auto subview = Kokkos::subview(view, Kokkos::ALL(), 0);
        for (std::size_t blockIdx = 0; blockIdx < x.size(); ++blockIdx)
        {
            if (!isOwned_[blockIdx]) continue;
            // Row map lists sub-rows of a block consecutively, so tRow0+j is valid
            const LocalIDType tRow0 = rowMap_->getLocalElement(
                static_cast<GlobalIDType>(localToGlobal_[blockIdx] * blockSize));
            for (std::size_t j = 0; j < blockSize; ++j)
                x[blockIdx][j] = subview(tRow0 + j);
        }
    }

    GridView gridView_;
    const DofMapper* mapper_ = nullptr;
    bool isParallel_;

    Teuchos::RCP<const Teuchos::Comm<int>> comm_;
    Teuchos::RCP<Teuchos::FancyOStream> fos_;
    Teuchos::RCP<const TpetraMap> rowMap_;
    Teuchos::RCP<const TpetraMap> colMap_;

    std::size_t rank_;
    std::size_t numBlocksGlobal_;
    std::vector<GlobalIDType> localToGlobal_;
    std::vector<bool> isOwned_;
    std::vector<bool> isGhost_;

    Teuchos::ParameterList params_;

    // Parallel MultiType: scalar-level and block-level DOF counts for each subdomain
    std::size_t nMomLocal_ = 0;       // scalar rows for momentum  (= momScalarBS_ * nMomBlockLocal_)
    std::size_t nMassLocal_ = 0;      // scalar rows for mass       (= massScalarBS_ * nMassBlockLocal_)
    std::size_t nMomBlockLocal_ = 0;  // block DOFs for momentum
    std::size_t nMassBlockLocal_ = 0; // block DOFs for mass
    // True when mass DOF rows are complete on the owner (FCDiamond-like, face-based).
    // False when mass rows need cross-process global assembly (Box-like, vertex-based).
    bool massRowsAreComplete_ = true;
    // Block-level closures (set by setupParallelMultiTypeFuncs_)
    std::function<void(MomBlockVec_&)>           momRhsPreprocess_;
    std::function<void(MassBlockVec_&)>          massRhsPreprocess_;
    std::function<void(MomBlockVec_&)>           momGhostBroadcastBlock_;
    std::function<void(MassBlockVec_&)>          massGhostBroadcastBlock_;
    std::function<void(std::vector<std::size_t>&)> momExtraCountComm_;
    std::function<void(std::vector<std::size_t>&)> massExtraCountComm_;
};

} // end namespace Dumux

#endif // DUMUX_HAVE_TRILINOS

#endif
