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

#include <dumux/assembly/jacobianpattern.hh>
#include <dumux/linear/parallelhelpers.hh>
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

/*!
 * \ingroup Linear
 * \brief Direct linear solvers from Trilinos Amesos2
 *
 * Wraps MUMPS (or other Amesos2 backends) via Tpetra for parallel direct solves.
 * Handles both sequential and parallel (non-overlapping and overlapping) decompositions.
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

    using GridView = typename LSTraits::GridView;
    using DofMapper = typename LSTraits::DofMapper;
    static constexpr int dofCodim = LSTraits::dofCodim;
    // Scalar DOFs per block DOF (1 for scalar PDEs, N for N-component systems)
    static constexpr std::size_t blockSize = XVector::block_type::size();

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

    template<class GridGeometry>
    DirectSolverAmesos2(const GridGeometry& gridGeometry,
                        const GridView& gridView,
                        const DofMapper& dofMapper,
                        const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
    , gridView_(gridView)
    , mapper_(dofMapper)
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

    SolverResult solve(const Matrix& A, XVector& x, const BVector& b)
    { return solve_(A, x, b); }

    Scalar norm(const XVector&) const { return 0.0; }

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
                VectorCommDataHandleMin<DofMapper, std::vector<GlobalIDType>, dofCodim, GlobalIDType>
                    handle(mapper, localToGlobal_);
                gridView.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);

                return;
            }
        }

        // Sequential: all DOFs are owned with consecutive IDs
        numBlocksGlobal_ = nLocal;
        std::fill(isOwned_.begin(), isOwned_.end(), true);
        for (std::size_t i = 0; i < nLocal; ++i)
            localToGlobal_[i] = static_cast<GlobalIDType>(i);
    }

    /*!
     * Row map: owned scalar DOFs only.
     * Scalar GID for (block i, sub k) = localToGlobal_[i] * blockSize + k.
     */
    Teuchos::RCP<const TpetraMap> buildRowMap_()
    {
        std::vector<GlobalIDType> ownedGIDs;
        ownedGIDs.reserve(
            static_cast<std::size_t>(std::count(isOwned_.begin(), isOwned_.end(), true)) * blockSize);
        for (std::size_t i = 0; i < mapper_.size(); ++i)
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
        const std::size_t nLocal = mapper_.size();
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
     * Uses the colMap invariant: Tpetra local col = blockCol*blockSize + subCol.
     */
    Teuchos::RCP<const Graph> buildGraph_(const Matrix& A)
    {
        const std::size_t nLocalTpetraRows = rowMap_->getLocalNumElements();

        // Count scalar entries per Tpetra row (each block row has rowIt->size() block cols,
        // each block col expands to blockSize scalar cols)
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

        auto graph = Teuchos::rcp(new Graph(rowMap_, colMap_, numEntPerRow()));

        // Insert column indices (using colMap invariant: local col = blockCol*BS + subCol)
        Teuchos::Array<LocalIDType> colIndices;
        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            const std::size_t blockRow = rowIt.index();
            if (!isOwned_[blockRow]) continue;

            colIndices.clear();
            colIndices.reserve(rowIt->size() * blockSize);
            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                for (std::size_t kj = 0; kj < blockSize; ++kj)
                    colIndices.push_back(
                        static_cast<LocalIDType>(colIt.index() * blockSize + kj));

            for (std::size_t ki = 0; ki < blockSize; ++ki)
            {
                const LocalIDType tRow = rowMap_->getLocalElement(
                    static_cast<GlobalIDType>(localToGlobal_[blockRow] * blockSize + ki));
                graph->insertLocalIndices(tRow, colIndices());
            }
        }

        graph->fillComplete(rowMap_, rowMap_);
        return graph;
    }

    SolverResult solve_(const Matrix& A, XVector& x, const BVector& b)
    {
        auto Acopy = A;
        auto bcopy = b;

        if constexpr (LSTraits::canCommunicate)
        {
            if (isParallel_)
            {
                if (LSTraits::isNonOverlapping(gridView_))
                {
                    // Non-overlapping decomposition (e.g. BOX):
                    // extend matrix pattern at border DOFs and accumulate cross-process contributions
                    using MatHelper = ParallelMatrixHelper<Matrix, GridView, DofMapper, dofCodim>;
                    MatHelper matrixHelper(gridView_, mapper_);
                    matrixHelper.extendMatrix(Acopy, [this](std::size_t idx){ return isGhost_[idx]; });
                    matrixHelper.sumEntries(Acopy);
                    ParallelVectorHelper<GridView, DofMapper, dofCodim> vecHelper(gridView_, mapper_);
                    vecHelper.makeNonOverlappingConsistent(bcopy);
                }
                // For overlapping decompositions (e.g. BOX with YaspGrid overlap=1):
                // The overlap layer already provides each process with a complete local
                // neighbourhood for all interior/border DOFs — no preprocessing needed.
                // Non-owned (overlap/ghost) DOF rows are simply skipped during Tpetra conversion.
            }
        }

        auto graph = buildGraph_(Acopy);
        auto AA = convertMatrix_(Acopy, graph);
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
                VectorCommDataHandleSum<DofMapper, XVector, dofCodim> handle(mapper_, x);
                gridView_.communicate(handle, Dune::All_All_Interface, Dune::ForwardCommunication);
            }
        }

        SolverResult result;
        result.converged = true;
        return result;
    }

    Teuchos::RCP<SolverMatrix> convertMatrix_(const Matrix& A, const Teuchos::RCP<const Graph>& graph)
    {
        Teuchos::RCP<SolverMatrix> AA = Teuchos::rcp(new SolverMatrix(graph));

        using BlockType = typename Matrix::block_type;
        static constexpr int BS_I = BlockType::rows;
        static constexpr int BS_J = BlockType::cols;

        Teuchos::Array<LocalIDType> cols(0);
        Teuchos::Array<Scalar> vals(0);

        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            const std::size_t blockRow = rowIt.index();
            if (!isOwned_[blockRow]) continue;

            const std::size_t nBlockCols = rowIt->size();
            cols.resize(nBlockCols * BS_J);
            vals.resize(nBlockCols * BS_J);

            for (int ki = 0; ki < BS_I; ++ki)
            {
                const LocalIDType tRow = rowMap_->getLocalElement(
                    static_cast<GlobalIDType>(localToGlobal_[blockRow] * BS_I + ki));

                int idx = 0;
                for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                {
                    const std::size_t blockCol = colIt.index();
                    for (int kj = 0; kj < BS_J; ++kj)
                    {
                        // colMap invariant: local col index = blockCol*BS_J + kj
                        cols[idx] = static_cast<LocalIDType>(blockCol * BS_J + kj);
                        vals[idx] = (*colIt)[ki][kj];
                        ++idx;
                    }
                }
                AA->replaceLocalValues(tRow, cols(), vals());
            }
        }

        AA->fillComplete();
        return AA;
    }

    Teuchos::RCP<SolverVector> convertBVector_(const BVector& b)
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

    void retrieveXVector_(const Teuchos::RCP<SolverVector>& XX, XVector& x)
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

    const GridView gridView_;
    const DofMapper& mapper_;
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
};

} // end namespace Dumux

#endif // DUMUX_HAVE_TRILINOS

#endif
