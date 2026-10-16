// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Fine-level smoothers of the p-multigrid preconditioner: an exact solve and restricted additive Schwarz
 *        smoothers on the subdomains of a parallel decomposition
 */
#ifndef DUMUX_LINEAR_PMULTIGRID_SMOOTHERS_HH
#define DUMUX_LINEAR_PMULTIGRID_SMOOTHERS_HH

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <memory>
#include <utility>
#include <vector>

#include <dune/common/fmatrixev.hh>
#include <dune/common/fvector.hh>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/solver.hh>
#include <dune/istl/solvercategory.hh>
#if HAVE_SUITESPARSE_UMFPACK
#include <dune/istl/umfpack.hh>
#endif
#if HAVE_MPI
#include <dune/istl/owneroverlapcopy.hh>
#endif

#include <dumux/common/multimapperview.hh>

namespace Dumux::Detail::PMultigrid {

#if HAVE_SUITESPARSE_UMFPACK
//! An exact solve with the LU factorization of the matrix, as smoother on the subdomain of a process
template<class M, class X, class Y>
class DirectSolve : public Dune::Preconditioner<X, Y>
{
public:
    explicit DirectSolve(const M& A) : solver_(A, 0) {}

    void pre(X&, Y&) override {}

    void apply(X& v, const Y& d) override
    {
        // the solver may overwrite the right-hand side
        rhs_ = d;
        Dune::InverseOperatorResult result;
        solver_.apply(v, rhs_, result);
    }

    void post(X&) override {}

    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::sequential; }

private:
    Dune::UMFPack<M> solver_;
    Y rhs_;
};
#endif

/*!
 * \brief Shift every singular value of a square block by m: D + m Q, with Q = U V^T the orthogonal factor of the
 *        singular value decomposition D = U S V^T
 * \note On the null space of a singular D, Q maps the right singular vectors to left singular vectors orthogonal to
 *       the range of D, which are the right singular vectors themselves where possible (Q is the identity on the null
 *       space of a symmetric D).
 * \note Blocks with more than three rows need LAPACK (Dune::FMatrixHelp::eigenValuesVectors).
 */
template<class Block, class Scalar>
void shiftSingularValues(Block& D, const Scalar m)
{
    constexpr int n = Block::rows;
    using Vector = Dune::FieldVector<Scalar, n>;
    Block DtD(0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k)
                DtD[i][j] += D[k][i]*D[k][j];

    Vector eigenValues;
    Block eigenVectors;
    Dune::FMatrixHelp::eigenValuesVectors(DtD, eigenValues, eigenVectors);
    const Scalar largest = eigenValues[n-1];

    // left singular vectors by decreasing singular value, orthonormalized (modified Gram-Schmidt)
    std::array<Vector, n> left;
    const auto orthogonalize = [&](Vector u, int numLeft)
    {
        for (int p = 0; p < numLeft; ++p)
            u.axpy(-(left[p]*u), left[p]);
        return u;
    };

    Block Q(0.0);
    for (int l = 0; l < n; ++l)
    {
        const int k = n-1-l;
        const Vector& v = eigenVectors[k];
        Vector u = v;
        if (largest > 0.0 && eigenValues[k] > 1e-14*largest)
            D.mv(v, u);
        u /= u.two_norm();
        u = orthogonalize(u, l);

        // a candidate (nearly) in the span of the previous left singular vectors is replaced by the
        // unit vector with the largest part orthogonal to that span (at least 1/sqrt(n) of its length)
        if (u.two_norm() < 0.5)
        {
            Scalar largestNorm = -1.0;
            for (int i = 0; i < n; ++i)
            {
                Vector unit(0.0);
                unit[i] = 1.0;
                const auto candidate = orthogonalize(unit, l);
                if (candidate.two_norm() > largestNorm)
                {
                    largestNorm = candidate.two_norm();
                    u = candidate;
                }
            }
        }

        u /= u.two_norm();
        left[l] = u;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                Q[i][j] += u[i]*v[j];
    }

    D.axpy(m, Q);
}

//! The sums of the absolute values of the entries of a matrix row, one per row of its blocks
template<class Row>
Dune::FieldVector<typename Row::block_type::field_type, Row::block_type::rows> absoluteRowSum(const Row& row)
{
    using Block = typename Row::block_type;
    Dune::FieldVector<typename Block::field_type, Block::rows> sum(0.0);
    for (auto col = row.begin(); col != row.end(); ++col)
        for (int r = 0; r < Block::rows; ++r)
            for (int c = 0; c < Block::cols; ++c)
                sum[r] += std::abs((*col)[r][c]);
    return sum;
}

/*!
 * \brief Enlarge the singular values of a diagonal block by the absolute values of the couplings its row lacks
 * \param D the diagonal block
 * \param completeRowSum the absolute row sums of the complete row
 * \param localRowSum the absolute row sums of the row as present
 * \note A diagonal block that stays zero becomes the identity.
 */
template<class Block, class RowSum>
void enlargeDiagonalBlock(Block& D, const RowSum& completeRowSum, const RowSum& localRowSum)
{
    typename Block::field_type missing = 0.0;
    for (int r = 0; r < Block::rows; ++r)
        missing = std::max(missing, completeRowSum[r] - localRowSum[r]);
    shiftSingularValues(D, missing);
    if (D.frobenius_norm() == 0.0)
        for (int k = 0; k < Block::rows; ++k)
            D[k][k] = 1.0;
}

/*!
 * \brief Enlarge the diagonal blocks of the rows marked in isIncomplete (see enlargeDiagonalBlock()), which keeps
 *        a subdomain problem with rows truncated at the subdomain boundary well posed
 * \param A the matrix of the subdomain
 * \param completeRowSums the absolute row sums of the complete rows, for the rows marked in isIncomplete
 * \param isIncomplete marks the rows that may lack couplings to dofs outside the subdomain
 */
template<class Matrix, class RowSums>
void enlargeIncompleteRows(Matrix& A, const RowSums& completeRowSums, const std::vector<bool>& isIncomplete)
{
    for (std::size_t i = 0; i < isIncomplete.size(); ++i)
        if (isIncomplete[i])
            enlargeDiagonalBlock(A[i][i], completeRowSums[i], absoluteRowSum(A[i]));
}

#if HAVE_MPI
//! Marks the dofs that the process does not own in the index set of a parallel communication
template<class Comm>
std::vector<bool> nonOwnedDofs(const Comm& communication, std::size_t numDofs)
{
    std::vector<bool> isNonOwned(numDofs, false);
    for (const auto& index : communication.indexSet())
        if (index.local().attribute() != Dune::OwnerOverlapCopyAttributeSet::owner)
            isNonOwned[index.local().local()] = true;
    return isNonOwned;
}
#endif

/*!
 * \brief Parallel smoother on the subdomain of a process including the dofs it does not own (restricted additive Schwarz)
 *
 * The local smoother acts on a matrix with the rows of all dofs of the process: overlap dofs of an overlapping
 * decomposition, or, for a nonoverlapping decomposition, ghost dofs whose rows are those of the owning processes.
 * The defect at the dofs not owned by the process is taken from the owners, and so is their correction.
 */
template<class X, class Y, class Comm>
class RestrictedSchwarzSmoother : public Dune::Preconditioner<X, Y>
{
public:
    RestrictedSchwarzSmoother(std::shared_ptr<Dune::Preconditioner<X, Y>> localSmoother,
                              std::shared_ptr<const void> localMatrix,
                              const Comm& communication,
                              Dune::SolverCategory::Category category)
    : localMatrix_(std::move(localMatrix))
    , localSmoother_(std::move(localSmoother))
    , communication_(communication)
    , category_(category)
    {}

    void pre(X& x, Y& b) override
    { localSmoother_->pre(x, b); }

    void apply(X& v, const Y& d) override
    {
        defect_ = d;
        communication_.copyOwnerToAll(defect_, defect_);
        localSmoother_->apply(v, defect_);
        communication_.copyOwnerToAll(v, v);
    }

    void post(X& x) override
    { localSmoother_->post(x); }

    Dune::SolverCategory::Category category() const override
    { return category_; }

private:
    // declared first, destroyed last: the local smoother refers to the local matrix
    std::shared_ptr<const void> localMatrix_;
    std::shared_ptr<Dune::Preconditioner<X, Y>> localSmoother_;
    const Comm& communication_;
    Dune::SolverCategory::Category category_;
    Y defect_;
};

//! Identifies a dof across processes: the entity carrying it and the position among the dofs of the entity
template<class IdType>
struct DofKey
{
    int codim;
    IdType id;
    int subIndex;

    friend bool operator<(const DofKey& a, const DofKey& b)
    {
        if (a.codim != b.codim)
            return a.codim < b.codim;
        if (a.id != b.id)
            return a.id < b.id;
        return a.subIndex < b.subIndex;
    }
};

/*!
 * \brief Extends the subdomain of a process by the dofs of all elements containing one of its border vertices
 *
 * A nonoverlapping decomposition provides at most the face neighbours of the elements of a process as ghost
 * elements, so the rows of border dofs lack the couplings to the dofs of elements that touch the subdomain at a
 * single vertex (or edge). At every border vertex, each process sends the rows of the dofs of its interior elements
 * containing the vertex. Received dofs not present on the receiving process get extended local indices
 * (after the local ones), their rows and the missing couplings of the local rows form the extended matrix.
 * The defect at the extended dofs is exchanged the same way. Dofs are identified by DofKey.
 */
template<class GridView, class Matrix, class Vector>
class VertexPatchOverlap
{
    using IdType = typename GridView::Grid::GlobalIdSet::IdType;
    using Key = DofKey<IdType>;
    using Block = typename Matrix::block_type;
    using VectorBlock = typename Vector::block_type;
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    static constexpr int dim = GridView::dimension;

    using Scalar = typename Block::field_type;
    using RowSums = Dune::FieldVector<Scalar, Block::rows>;
    struct MatrixItem { Key row; Key col; Block value; };
    struct VectorItem { Key dof; VectorBlock value; };
    struct RowSumItem { Key dof; RowSums value; };

    //! Sends the items of the dofs around the border vertices, scattered into a callback
    template<class Item, class Gather, class Scatter>
    class VertexExchange : public Dune::CommDataHandleIF<VertexExchange<Item, Gather, Scatter>, Item>
    {
    public:
        VertexExchange(const VertexMapper& vertexMapper, Gather gather, Scatter scatter, std::vector<std::size_t> numItems)
        : vertexMapper_(vertexMapper), gather_(std::move(gather)), scatter_(std::move(scatter)), numItems_(std::move(numItems))
        {}

        bool contains(int, int codim) const { return codim == dim; }
        bool fixedSize(int, int) const { return false; }

        template<class Entity>
        std::size_t size(const Entity& e) const { return numItems_[vertexMapper_.index(e)]; }

        template<class Buffer, class Entity>
        void gather(Buffer& buffer, const Entity& e) const
        { gather_(vertexMapper_.index(e), [&](const Item& item){ buffer.write(item); }); }

        template<class Buffer, class Entity>
        void scatter(Buffer& buffer, const Entity&, std::size_t n)
        {
            for (std::size_t k = 0; k < n; ++k)
            {
                Item item;
                buffer.read(item);
                scatter_(item);
            }
        }

    private:
        const VertexMapper& vertexMapper_;
        Gather gather_;
        Scatter scatter_;
        std::vector<std::size_t> numItems_;
    };

public:
    template<class GridGeometry>
    VertexPatchOverlap(const GridGeometry& gridGeometry, const VertexMapper& vertexMapper)
    : gridView_(gridGeometry.gridView())
    , vertexMapper_(vertexMapper)
    , numLocal_(gridGeometry.numDofs())
    , keys_(numLocal_)
    , patchDofs_(vertexMapper.size())
    {
        const auto& idSet = gridView_.grid().globalIdSet();
        const auto& dofMapper = gridGeometry.dofMapper();
        for (const auto& element : elements(gridView_))
        {
            Dune::Hybrid::forEach(std::make_integer_sequence<int, dim+1>{}, [&](auto codim)
            {
                for (int i = 0; i < static_cast<int>(element.subEntities(codim)); ++i)
                {
                    const auto entity = element.template subEntity<codim>(i);
                    int subIndex = 0;
                    for (const auto dofIdx : asMultiMapper(dofMapper).indices(entity))
                        keys_[dofIdx] = Key{static_cast<int>(codim), idSet.id(entity), subIndex++};
                }
            });
        }
        for (std::size_t i = 0; i < numLocal_; ++i)
            localIndex_.emplace(keys_[i], i);

        auto elemDisc = localView(gridGeometry);
        for (const auto& element : elements(gridView_, Dune::Partitions::interior))
        {
            bool bound = false;
            for (int i = 0; i < static_cast<int>(element.subEntities(dim)); ++i)
            {
                const auto vertex = element.template subEntity<dim>(i);
                if (vertex.partitionType() != Dune::BorderEntity)
                    continue;
                if (!bound)
                {
                    elemDisc.bind(element);
                    bound = true;
                }
                auto& patch = patchDofs_[vertexMapper_.index(vertex)];
                for (const auto& localDof : localDofs(elemDisc))
                    patch.push_back(localDof.dofIndex());
            }
        }
        for (auto& patch : patchDofs_)
        {
            std::sort(patch.begin(), patch.end());
            patch.erase(std::unique(patch.begin(), patch.end()), patch.end());
        }
    }

    // referenced by the communication of every application
    VertexPatchOverlap(const VertexPatchOverlap&) = delete;
    VertexPatchOverlap& operator=(const VertexPatchOverlap&) = delete;

    //! The subdomain extended for one matrix: its matrix and the indices of the extended dofs (after the local ones)
    struct Extension
    {
        Matrix matrix;
        std::map<Key, std::size_t> index;
    };

    /*!
     * \brief The matrix with the rows and columns of the local and the extended dofs (collective)
     * \param localMatrix the local matrix, whose entries between local dofs are complete
     * \param extension the extended subdomain (output), with a matrix not yet built
     * \note The rows of the extended dofs lack the couplings to dofs outside the extended subdomain; the absolute
     *       values of these couplings enlarge the singular values of their diagonal blocks (enlargeDiagonalBlock()),
     *       which keeps the extended subdomain problem well posed.
     */
    void extend(const Matrix& localMatrix, Extension& extension) const
    {
        auto& extendedIndex = extension.index;
        std::vector<std::size_t> numItems(patchDofs_.size(), 0);
        for (std::size_t v = 0; v < patchDofs_.size(); ++v)
            for (const auto dofIdx : patchDofs_[v])
                numItems[v] += localMatrix[dofIdx].size();

        std::vector<MatrixItem> received;
        auto gather = [&](std::size_t v, const auto& write)
        {
            for (const auto dofIdx : patchDofs_[v])
                for (auto col = localMatrix[dofIdx].begin(); col != localMatrix[dofIdx].end(); ++col)
                    write(MatrixItem{keys_[dofIdx], keys_[col.index()], *col});
        };
        auto scatter = [&](const MatrixItem& item){ received.push_back(item); };
        VertexExchange<MatrixItem, decltype(gather), decltype(scatter)> exchange(vertexMapper_, gather, scatter, std::move(numItems));
        gridView_.communicate(exchange, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);

        // the received rows of dofs not present on this process get extended indices
        extendedIndex.clear();
        for (const auto& item : received)
            if (!localIndex_.count(item.row) && !extendedIndex.count(item.row))
                extendedIndex.emplace(item.row, numLocal_ + extendedIndex.size());
        const std::size_t numDofs = numLocal_ + extendedIndex.size();

        const auto index = [&](const Key& key) -> std::ptrdiff_t
        {
            if (const auto it = localIndex_.find(key); it != localIndex_.end())
                return it->second;
            if (const auto it = extendedIndex.find(key); it != extendedIndex.end())
                return it->second;
            return -1;
        };

        // entries between local dofs are complete in the local matrix; the others come from the received rows
        std::map<std::pair<std::size_t, std::size_t>, Block> added;
        for (const auto& item : received)
        {
            const auto row = index(item.row);
            const auto col = index(item.col);
            if (row < 0 || col < 0 || (static_cast<std::size_t>(row) < numLocal_ && static_cast<std::size_t>(col) < numLocal_))
                continue;
            added[{row, col}] = item.value;
        }

        std::vector<std::vector<std::size_t>> pattern(numDofs);
        for (auto row = localMatrix.begin(); row != localMatrix.end(); ++row)
            for (auto col = row->begin(); col != row->end(); ++col)
                pattern[row.index()].push_back(col.index());
        for (const auto& [entry, value] : added)
            pattern[entry.first].push_back(entry.second);
        for (std::size_t i = numLocal_; i < numDofs; ++i)
            pattern[i].push_back(i);
        std::size_t numNonZeros = 0;
        for (auto& cols : pattern)
        {
            std::sort(cols.begin(), cols.end());
            cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
            numNonZeros += cols.size();
        }

        auto& extended = extension.matrix;
        extended.setBuildMode(Matrix::row_wise);
        extended.setSize(numDofs, numDofs, numNonZeros);
        for (auto row = extended.createbegin(); row != extended.createend(); ++row)
            for (const auto col : pattern[row.index()])
                row.insert(col);
        extended = 0.0;
        for (auto row = localMatrix.begin(); row != localMatrix.end(); ++row)
            for (auto col = row->begin(); col != row->end(); ++col)
                extended[row.index()][col.index()] = *col;
        for (const auto& [entry, value] : added)
            extended[entry.first][entry.second] = value;
        // the complete absolute row sums of the extended dofs, from the processes sending their rows
        std::map<Key, RowSums> completeRowSums;
        {
            std::vector<std::size_t> numSums(patchDofs_.size());
            for (std::size_t v = 0; v < patchDofs_.size(); ++v)
                numSums[v] = patchDofs_[v].size();
            auto gatherSums = [&](std::size_t v, const auto& write)
            {
                for (const auto dofIdx : patchDofs_[v])
                    write(RowSumItem{keys_[dofIdx], absoluteRowSum(localMatrix[dofIdx])});
            };
            auto scatterSums = [&](const RowSumItem& item)
            {
                if (extendedIndex.count(item.dof))
                {
                    auto& sums = completeRowSums[item.dof];
                    for (int r = 0; r < Block::rows; ++r)
                        sums[r] = std::max(sums[r], item.value[r]);
                }
            };
            VertexExchange<RowSumItem, decltype(gatherSums), decltype(scatterSums)> sumExchange(
                vertexMapper_, gatherSums, scatterSums, std::move(numSums)
            );
            gridView_.communicate(sumExchange, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
        }

        for (const auto& [key, i] : extendedIndex)
        {
            const auto localSums = absoluteRowSum(extended[i]);
            const auto it = completeRowSums.find(key);
            const auto& completeSums = it != completeRowSums.end() ? it->second : localSums;
            enlargeDiagonalBlock(extended[i][i], completeSums, localSums);
        }
    }

    /*!
     * \brief The defect at the local and the extended dofs (collective)
     * \param extension the extended subdomain
     * \param localDefect the defect at the local dofs, consistent
     * \param extendedDefect the defect at the local and extended dofs (output)
     */
    void extendDefect(const Extension& extension, const Vector& localDefect, Vector& extendedDefect) const
    {
        const auto& extendedIndex = extension.index;
        extendedDefect.resize(numLocal_ + extendedIndex.size());
        extendedDefect = 0.0;
        for (std::size_t i = 0; i < numLocal_; ++i)
            extendedDefect[i] = localDefect[i];

        std::vector<std::size_t> numItems(patchDofs_.size());
        for (std::size_t v = 0; v < patchDofs_.size(); ++v)
            numItems[v] = patchDofs_[v].size();
        auto gather = [&](std::size_t v, const auto& write)
        {
            for (const auto dofIdx : patchDofs_[v])
                write(VectorItem{keys_[dofIdx], localDefect[dofIdx]});
        };
        auto scatter = [&](const VectorItem& item)
        {
            if (const auto it = extendedIndex.find(item.dof); it != extendedIndex.end())
                extendedDefect[it->second] = item.value;
        };
        VertexExchange<VectorItem, decltype(gather), decltype(scatter)> exchange(vertexMapper_, gather, scatter, std::move(numItems));
        gridView_.communicate(exchange, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
    }

private:
    GridView gridView_;
    const VertexMapper& vertexMapper_;
    std::size_t numLocal_;
    std::vector<Key> keys_;
    std::map<Key, std::size_t> localIndex_;
    std::vector<std::vector<std::size_t>> patchDofs_;
};

/*!
 * \brief Parallel smoother on the subdomain of a process extended by the vertex patches of its border dofs
 *        (restricted additive Schwarz)
 *
 * The local smoother acts on the matrix of a VertexPatchOverlap::Extension. The defect at the ghost dofs and at the
 * extended dofs is taken from the processes owning them, and so is the correction of the dofs not owned by the process.
 */
template<class X, class Y, class Comm, class Overlap>
class VertexPatchSchwarzSmoother : public Dune::Preconditioner<X, Y>
{
    using Extension = typename Overlap::Extension;
public:
    VertexPatchSchwarzSmoother(std::shared_ptr<Dune::Preconditioner<X, Y>> localSmoother,
                               std::shared_ptr<const Extension> extension,
                               std::shared_ptr<const Overlap> overlap,
                               const Comm& communication)
    : extension_(std::move(extension))
    , localSmoother_(std::move(localSmoother))
    , overlap_(std::move(overlap))
    , communication_(communication)
    {}

    void pre(X&, Y&) override
    {
        extendedCorrection_.resize(extension_->matrix.N());
        extendedCorrection_ = 0.0;
        extendedDefect_.resize(extension_->matrix.N());
        extendedDefect_ = 0.0;
        localSmoother_->pre(extendedCorrection_, extendedDefect_);
    }

    void apply(X& v, const Y& d) override
    {
        defect_ = d;
        communication_.copyOwnerToAll(defect_, defect_);
        overlap_->extendDefect(*extension_, defect_, extendedDefect_);
        extendedCorrection_.resize(extendedDefect_.size());
        extendedCorrection_ = 0.0;
        localSmoother_->apply(extendedCorrection_, extendedDefect_);
        for (std::size_t i = 0; i < v.size(); ++i)
            v[i] = extendedCorrection_[i];
        communication_.copyOwnerToAll(v, v);
    }

    void post(X&) override
    { localSmoother_->post(extendedCorrection_); }

    Dune::SolverCategory::Category category() const override
    { return Dune::SolverCategory::nonoverlapping; }

private:
    // declared first, destroyed last: the local smoother refers to the extended matrix
    std::shared_ptr<const Extension> extension_;
    std::shared_ptr<Dune::Preconditioner<X, Y>> localSmoother_;
    std::shared_ptr<const Overlap> overlap_;
    const Comm& communication_;
    Y defect_, extendedDefect_;
    X extendedCorrection_;
};

} // end namespace Dumux::Detail::PMultigrid

#endif
