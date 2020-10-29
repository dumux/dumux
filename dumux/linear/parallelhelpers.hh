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
 * \brief Provides a helper class for nonoverlapping
 *        decomposition.
 */
#ifndef DUMUX_LINEAR_PARALLELHELPERS_HH
#define DUMUX_LINEAR_PARALLELHELPERS_HH

#if HAVE_MPI

#include <dune/geometry/dimension.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A parallel helper class providing a nonoverlapping
 *        decomposition of all degrees of freedom
 */
// operator that resets result to zero at constrained DOFS
template<class LinearSolverTraits>
class ParallelISTLHelper
{
    using GridView = typename LinearSolverTraits::GridView;
    using DofMapper = typename LinearSolverTraits::DofMapper;
    enum { dofCodim = LinearSolverTraits::dofCodim };

    // TODO: this is some large number (replace by limits?)
    static constexpr std::size_t ghostMarker_ = 1<<24;

    class BaseGatherScatter
    {
    public:
        BaseGatherScatter(const DofMapper& mapper)
        : mapper_(mapper) {}

        template<class EntityType>
        int index(const EntityType& e) const
        { return mapper_.index(e); }

        bool contains(int dim, int codim) const
        { return dofCodim == codim; }

        //! returns true if size per entity of given dim and codim is a constant
        bool fixedSize(int dim, int codim) const
        { return true; }

        template<class EntityType>
        std::size_t size(EntityType& e) const
        { return 1; }

        template<class EntityType>
        bool isNeitherInteriorNorBorderEntity(EntityType& e) const
        { return e.partitionType() != Dune::InteriorEntity && e.partitionType() != Dune::BorderEntity; }

    private:
        const DofMapper& mapper_;
    };

    /*!
     * \brief Writes ghostMarker_ to each data item (of the container) that is gathered or scattered
     * and is neither interior nor border.
     *
     * Can be used to mark ghost cells.
     */
    class GhostGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<GhostGatherScatter,std::size_t>
    {
    public:
        using DataType = std::size_t;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedSize;
        using BaseGatherScatter::size;

        GhostGatherScatter(std::vector<std::size_t>& ranks, const DofMapper& mapper)
        : BaseGatherScatter(mapper), ranks_(ranks)
        {}

        template<class MessageBuffer, class EntityType>
        void gather(MessageBuffer& buff, const EntityType& e) const
        {
            auto& data = ranks_[this->index(e)];
            if (this->isNeitherInteriorNorBorderEntity(e))
                data = ghostMarker_;
            buff.write(data);
        }

        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, std::size_t n)
        {
            std::size_t x;
            buff.read(x);
            auto& data = ranks_[this->index(e)];
            if (this->isNeitherInteriorNorBorderEntity(e))
                data = ghostMarker_;
        }
    private:
        std::vector<std::size_t>& ranks_;
    };

    /*!
     * \brief GatherScatter handle that sets ghostMarker_ for data items neither associated to
     * the interior or border and take the minimum when scattering.
     *
     * Used to compute an owner rank for each unknown.
     */
    class InteriorBorderGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<InteriorBorderGatherScatter,std::size_t>
    {
    public:
        using DataType = std::size_t;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedSize;
        using BaseGatherScatter::size;

        InteriorBorderGatherScatter(std::vector<std::size_t>& ranks, const DofMapper& mapper)
        : BaseGatherScatter(mapper), ranks_(ranks)
        {}

        template<class MessageBuffer, class EntityType>
        void gather(MessageBuffer& buff, const EntityType& e) const
        {
            auto& data = ranks_[this->index(e)];
            if (this->isNeitherInteriorNorBorderEntity(e))
                data = ghostMarker_;
            buff.write(data);
        }

        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, std::size_t n)
        {
            std::size_t x;
            buff.read(x);
            auto& data = ranks_[this->index(e)];
            using std::min;
            data = this->isNeitherInteriorNorBorderEntity(e) ? x : min(data,x);
        }
    private:
        std::vector<std::size_t>& ranks_;
    };

    /*!
     * \brief GatherScatter handle for finding out about neighbouring processor ranks.
     *
     */
    struct NeighbourGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<NeighbourGatherScatter,int>
    {
        using DataType = int;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedSize;
        using BaseGatherScatter::size;

        NeighbourGatherScatter(const DofMapper& mapper, int rank, std::set<int>& neighbours)
        : BaseGatherScatter(mapper), rank_(rank), neighbours_(neighbours)
        {}

        template<class MessageBuffer, class EntityType>
        void gather(MessageBuffer& buff, const EntityType& e) const
        {
            buff.write(rank_);
        }

        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, std::size_t n)
        {
            int x;
            buff.read(x);
            neighbours_.insert(x);
        }
    private:
        int rank_;
        std::set<int>& neighbours_;
    };


    /*!
     * \brief GatherScatter handle for finding out about neighbouring processor ranks.
     *
     */
    struct SharedGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<SharedGatherScatter,int>
    {
        using DataType = int;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedSize;
        using BaseGatherScatter::size;

        SharedGatherScatter(std::vector<int>& shared, const DofMapper& mapper)
        : BaseGatherScatter(mapper), shared_(shared)
        {}

        template<class MessageBuffer, class EntityType>
        void gather(MessageBuffer& buff, EntityType& e) const
        {
            int data = true;
            buff.write(data);
        }

        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType &e, std::size_t n)
        {
            int x;
            buff.read(x);
            auto& data = shared_[this->index(e)];
            data = data || x;
        }
    private:
        std::vector<int>& shared_;
    };

    /*!
     * \brief GatherScatter handle for finding out about neighbouring processor ranks.
     *
     */
    template<typename GlobalIndex>
    struct GlobalIndexGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<GlobalIndexGatherScatter<GlobalIndex>, GlobalIndex>
    {
        using DataType = GlobalIndex;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedSize;
        using BaseGatherScatter::size;

        GlobalIndexGatherScatter(std::vector<GlobalIndex>& globalIndices, const DofMapper& mapper)
        : BaseGatherScatter(mapper), globalIndices_(globalIndices)
        {}

        template<class MessageBuffer, class EntityType>
        void gather(MessageBuffer& buff, const EntityType& e) const
        {
            buff.write(globalIndices_[this->index(e)]);
        }

        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, std::size_t n)
        {
            DataType x;
            buff.read(x);
            using std::min;
            auto& data = globalIndices_[this->index(e)];
            data = min(data, x);
        }
    private:
        std::vector<GlobalIndex>& globalIndices_;
    };

public:

    ParallelISTLHelper(const GridView& gridView, const DofMapper& mapper)
    : gridView_(gridView), mapper_(mapper)
    {
        initGhostsAndOwners_();
    }

    bool isGhost(std::size_t i) const
    { return isGhost_[i] == ghostMarker_; }

    /*!
     * \brief Creates a parallel index set
     *
     * \tparam Comm The type of the OwnerOverlapCopyCommunication
     * communicators.
     */
    template<class Comm>
    void createParallelIndexSet(Comm& comm) const
    {
        if (gridView_.comm().size() <= 1)
        {
            comm.remoteIndices().template rebuild<false>();
            return;
        }

        // First find out which dofs we share with other processors
        std::vector<int> isShared(mapper_.size(), false);

        SharedGatherScatter sgs(isShared, mapper_);
        gridView_.communicate(sgs, Dune::All_All_Interface, Dune::ForwardCommunication);

        // Count shared dofs that we own
        using GlobalIndex = typename Comm::ParallelIndexSet::GlobalIndex;
        GlobalIndex count = 0;
        for (std::size_t i = 0; i < isShared.size(); ++i)
            if (isShared[i] && isOwned_[i] == 1)
                ++count;

        std::vector<GlobalIndex> counts(gridView_.comm().size());
        gridView_.comm().allgather(&count, 1, counts.data());

        // Compute start index start_p = \sum_{i=0}^{i<p} counts_i
        const int rank = gridView_.comm().rank();
        auto start = std::accumulate(counts.begin(), counts.begin() + rank, GlobalIndex(0));

        std::vector<GlobalIndex> globalIndices(mapper_.size(), std::numeric_limits<GlobalIndex>::max());

        for (std::size_t i = 0; i < globalIndices.size(); ++i)
        {
            if (isOwned_[i] == 1 && isShared[i])
            {
                globalIndices[i] = start; // GlobalIndex does not implement postfix ++
                ++start;
            }
        }

        // publish global indices for the shared DOFS to other processors.
        GlobalIndexGatherScatter<GlobalIndex> gigs(globalIndices, mapper_);
        gridView_.communicate(gigs, Dune::All_All_Interface, Dune::ForwardCommunication);

        resizeIndexSet_(comm, globalIndices);

        // Compute neighbours using communication
        std::set<int> neighbours;
        NeighbourGatherScatter ngs(mapper_, gridView_.comm().rank(), neighbours);
        gridView_.communicate(ngs, Dune::All_All_Interface, Dune::ForwardCommunication);
        comm.remoteIndices().setNeighbours(neighbours);

        comm.remoteIndices().template rebuild<false>();
    }

    //! Return the dofMapper
    const DofMapper& dofMapper() const
    { return mapper_; }

    //! Return the gridView
    const GridView& gridView() const
    { return gridView_; }

private:
    void initGhostsAndOwners_()
    {
        const auto rank = gridView_.comm().rank();
        isOwned_.assign(mapper_.size(), rank);
        // find out about ghosts
        GhostGatherScatter ggs(isOwned_, mapper_);

        if (gridView_.comm().size() > 1)
            gridView_.communicate(ggs, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);

        isGhost_ = isOwned_;

        // partition interior/border
        InteriorBorderGatherScatter dh(isOwned_, mapper_);

        if (gridView_.comm().size() > 1)
            gridView_.communicate(dh, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);

        // convert vector into mask vector
        for (auto& v : isOwned_)
            v = (v == rank) ? 1 : 0;
    }

    template<class Comm, class GlobalIndices>
    void resizeIndexSet_(Comm& comm, const GlobalIndices& globalIndices) const
    {
        comm.indexSet().beginResize();

        for (std::size_t localIndex = 0; localIndex < globalIndices.size(); ++localIndex)
        {
            const auto globalIndex = globalIndices[localIndex];
            if (globalIndex != std::numeric_limits<typename GlobalIndices::value_type>::max())
            {
                const bool isOwned = isOwned_[localIndex] > 0;
                const auto attr = getAttribute_(comm, isOwned, isGhost(localIndex));
                using LocalIndex = typename Comm::ParallelIndexSet::LocalIndex;
                comm.indexSet().add(globalIndex, LocalIndex{localIndex, attr});
            }
        }

        comm.indexSet().endResize();
    }

    template<class Comm>
    Dune::OwnerOverlapCopyAttributeSet::AttributeSet
    getAttribute_(const Comm& comm, const bool isOwned, const bool isGhost) const
    {
        if (isOwned)
            return Dune::OwnerOverlapCopyAttributeSet::owner;
        else if (isGhost && (comm.category() == static_cast<int>(Dune::SolverCategory::nonoverlapping)) )
            return Dune::OwnerOverlapCopyAttributeSet::overlap;
        else
            return Dune::OwnerOverlapCopyAttributeSet::copy;
    }

    const GridView gridView_; //!< the grid view
    const DofMapper& mapper_; //!< the dof mapper
    std::vector<std::size_t> isOwned_; //!< vector to identify unique decomposition
    std::vector<std::size_t> isGhost_; //!< vector to identify ghost dofs

}; // class ParallelISTLHelper

template<class GridView, class DofMapper, int dofCodim>
class ParallelVectorHelper
{
public:
    ParallelVectorHelper(const GridView& gridView, const DofMapper& mapper)
    : gridView_(gridView), mapper_(mapper)
    {}

    // \brief Make a vector of the box model consistent.
    template<class Block, class Alloc>
    void makeNonOverlappingConsistent(Dune::BlockVector<Block, Alloc>& v) const
    {
        VectorCommDataHandleSum<DofMapper, Dune::BlockVector<Block, Alloc>, dofCodim, Block> gs(mapper_, v);
        if (gridView_.comm().size() > 1)
            gridView_.communicate(gs, Dune::InteriorBorder_InteriorBorder_Interface,
                                  Dune::ForwardCommunication);
    }

private:
    const GridView gridView_; //!< the grid view
    const DofMapper& mapper_; //!< the dof mapper
};

/*!
 * \ingroup Linear
 * \brief Helper class for adding up matrix entries on border.
 */
template<class Matrix, class GridView, class DofMapper, int dofCodim>
class ParallelMatrixHelper
{
    static constexpr int dim = GridView::dimension;
    using Grid = typename GridView::Traits::Grid;
    using BlockType = typename Matrix::block_type;
    using IDS = typename Grid::Traits::GlobalIdSet;
    using IdType = typename IDS::IdType;

    /*!
     * \brief A DataHandle class to exchange matrix sparsity patterns.
     *
     *  We look at a 2D example with a nonoverlapping grid,
     *  two processes and no ghosts with Q1 discretization.
     *  Process 0 has the left part of the domain
     *  with three cells and eight vertices (1-8),
     *  Process 1 the right part with three cells
     *  and eight vertices (2,4,7-12).
     *  <pre>
     *  1 _ 2        2 _ 9 _ 10
     *  |   |        |   |   |
     *  3 _ 4 _ 7    4 _ 7 _ 11
     *  |   |   |        |   |
     *  5 _ 6 _ 8        8 _ 12
     *  </pre>
     *  If we look at vertex 7 and the corresponding entries in the matrix for P0,
     *  there will be entries for (7,4) and (7,8), but not for (7,2).
     *  The MatrixPatternExchange class will find these entries and fills a vector "sparsity",
     *  that contains all missing connections.
     */
    template<class IsGhostFunc>
    struct MatrixPatternExchange
    : public Dune::CommDataHandleIF<MatrixPatternExchange<IsGhostFunc>, IdType>
    {
        //! Export type of data for message buffer
        using DataType = IdType;

        MatrixPatternExchange(const DofMapper& mapper,
                              const std::map<IdType,int>& globalToLocal,
                              const std::map<int,IdType>& localToGlobal,
                              Matrix& A,
                              std::vector<std::set<int>>& sparsity,
                              const IsGhostFunc& isGhost)
        : mapper_(mapper), idToIndex_(globalToLocal), indexToID_(localToGlobal)
        , sparsity_(sparsity), A_(A), isGhost_(isGhost)
        {
            sparsity_.resize(A.N());
        }

        /*!
         * \brief Returns true if data for given valid codim should be communicated
         */
        bool contains (int dim, int codim) const
        { return (codim == dofCodim); }

        //! returns true if size per entity of given dim and codim is a constant
        bool fixedSize(int dim, int codim) const
        { return false; }

        /*!
         * \brief How many objects of type DataType have to be sent for a given entity
         */
        template<class EntityType>
        std::size_t size(EntityType& e) const
        {
            const auto rowIdx = mapper_.index(e);
            std::size_t  n = 0;
            for (auto colIt = A_[rowIdx].begin(); colIt != A_[rowIdx].end(); ++colIt)
                if (indexToID_.count(colIt.index()))
                    n++;

            return n;
        }

        /*!
         * \brief Pack data from user to message buffer
         */
        template<class MessageBuffer, class EntityType>
        void gather(MessageBuffer& buff, const EntityType& e) const
        {
            const auto rowIdx = mapper_.index(e);
            for (auto colIt = A_[rowIdx].begin(); colIt != A_[rowIdx].end(); ++colIt)
            {
                auto it = indexToID_.find(colIt.index());
                if (it != indexToID_.end())
                    buff.write(it->second);
            }
        }

        /*!
         * \brief Unpack data from message buffer to user
         */
        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, std::size_t n)
        {
            const auto rowIdx = mapper_.index(e);
            for (std::size_t k = 0; k < n; k++)
            {
                IdType id;
                buff.read(id);
                // only add entries corresponding to border entities
                const auto it = idToIndex_.find(id);
                if (it != idToIndex_.end())
                {
                    const auto colIdx = it->second;
                    if (!sparsity_[rowIdx].count(colIdx) && !isGhost_(colIdx))
                        sparsity_[rowIdx].insert(colIdx);
                }
            }
        }

    private:
        const DofMapper& mapper_;
        const std::map<IdType,int>& idToIndex_;
        const std::map<int,IdType>& indexToID_;
        std::vector<std::set<int>>& sparsity_;
        Matrix& A_;
        const IsGhostFunc& isGhost_;

    }; // class MatrixPatternExchange

    //! Local matrix blocks associated with the global id set
    struct MatrixEntry
    {
        IdType first;
        BlockType second;
        MatrixEntry (const IdType& f, const BlockType& s) : first(f), second(s) {}
        MatrixEntry () {}
    };

    //! A DataHandle class to exchange matrix entries
    struct MatrixEntryExchange
    : public Dune::CommDataHandleIF<MatrixEntryExchange, MatrixEntry>
    {
        //! Export type of data for message buffer
        using DataType = MatrixEntry;

        MatrixEntryExchange(const DofMapper& mapper,
                            const std::map<IdType,int>& globalToLocal,
                            const std::map<int,IdType>& localToGlobal,
                            Matrix& A)
        : mapper_(mapper), idToIndex_(globalToLocal), indexToID_(localToGlobal), A_(A)
        {}

        /*!
         * \brief Returns true if data for given valid codim should be communicated
         */
        bool contains(int dim, int codim) const
        { return (codim == dofCodim); }

        //! returns true if size per entity of given dim and codim is a constant
        bool fixedSize(int dim, int codim) const
        { return false; }

        /*!
         * \brief How many objects of type DataType have to be sent for a given entity
         */
        template<class EntityType>
        std::size_t size(EntityType& e) const
        {
            const auto rowIdx = mapper_.index(e);
            std::size_t n = 0;
            for (auto colIt = A_[rowIdx].begin(); colIt != A_[rowIdx].end(); ++colIt)
                if (indexToID_.count(colIt.index()))
                    n++;

            return n;
        }

        /*!
         * \brief Pack data from user to message buffer
         */
        template<class MessageBuffer, class EntityType>
        void gather(MessageBuffer& buff, const EntityType& e) const
        {
            const auto rowIdx = mapper_.index(e);
            for (auto colIt = A_[rowIdx].begin(); colIt != A_[rowIdx].end(); ++colIt)
            {
                auto it = indexToID_.find(colIt.index());
                if (it != indexToID_.end())
                    buff.write(MatrixEntry(it->second,*colIt));
            }
        }

        /*!
         * \brief Unpack data from message buffer to user
         */
        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, std::size_t n)
        {
            const auto rowIdx = mapper_.index(e);
            for (std::size_t k = 0; k < n; k++)
            {
                MatrixEntry m;
                buff.read(m);
                // only add entries corresponding to border entities
                auto it = idToIndex_.find(m.first);
                if (it != idToIndex_.end())
                    if (A_[rowIdx].find(it->second) != A_[rowIdx].end())
                        A_[rowIdx][it->second] += m.second;
            }
        }

    private:
        const DofMapper& mapper_;
        const std::map<IdType,int>& idToIndex_;
        const std::map<int,IdType>& indexToID_;
        Matrix& A_;

    }; // class MatrixEntryExchange

public:

    ParallelMatrixHelper(const GridView& gridView, const DofMapper& mapper)
    : gridView_(gridView), mapper_(mapper)
    {
        idToIndex_.clear();
        indexToID_.clear();

        for (const auto& entity : entities(gridView_, Dune::Codim<dofCodim>()))
        {
            if (entity.partitionType() == Dune::BorderEntity)
            {
                const int localIdx = mapper_.index(entity);
                IdType dofIdxGlobal = gridView_.grid().globalIdSet().id(entity);

                idToIndex_.emplace(dofIdxGlobal, localIdx);
                indexToID_.emplace(localIdx, dofIdxGlobal);
            }
        }
    }

    /*!
     * \brief communicates values for the sparsity pattern of the new matrix.
     * \param A Matrix to operate on.
     * \param isGhost Function returning if is ghost.
     */
    template<class IsGhostFunc>
    void extendMatrix(Matrix& A, const IsGhostFunc& isGhost)
    {
        if (gridView_.comm().size() <= 1)
            return;

        Matrix tmp(A);
        std::size_t nnz = 0;
        // get entries from other processes
        std::vector<std::set<int>> sparsity;
        MatrixPatternExchange<IsGhostFunc> datahandle(mapper_, idToIndex_, indexToID_, A, sparsity, isGhost);
        gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface,
                                          Dune::ForwardCommunication);
        // add own entries, count number of nonzeros
        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            for (auto colIt = A[rowIt.index()].begin(); colIt != A[rowIt.index()].end(); ++colIt)
                if (!sparsity[rowIt.index()].count(colIt.index()))
                    sparsity[rowIt.index()].insert(colIt.index());

            nnz += sparsity[rowIt.index()].size();
        }

        A.setSize(tmp.N(), tmp.N(), nnz);
        A.setBuildMode(Matrix::row_wise);
        auto citer = A.createbegin();
        for (auto i = sparsity.begin(), end = sparsity.end(); i!=end; ++i, ++citer)
        {
            for (auto si = i->begin(), send = i->end(); si!=send; ++si)
                citer.insert(*si);
        }

        // set matrix old values
        A = 0;
        for (auto rowIt = tmp.begin(); rowIt != tmp.end(); ++rowIt)
            for (auto colIt = tmp[rowIt.index()].begin(); colIt != tmp[rowIt.index()].end(); ++colIt)
                A[rowIt.index()][colIt.index()] = tmp[rowIt.index()][colIt.index()];
    }

    /*!
     * \brief Sums up the entries corresponding to border vertices.
     * \param A Matrix to operate on.
     */
    void sumEntries(Matrix& A)
    {
        if (gridView_.comm().size() <= 1)
            return;

        MatrixEntryExchange datahandle(mapper_, idToIndex_, indexToID_, A);
        gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface,
                                          Dune::ForwardCommunication);
    }

private:
    const GridView gridView_;
    const DofMapper& mapper_;
    std::map<IdType, int> idToIndex_;
    std::map<int, IdType> indexToID_;

};

/*!
 * \brief Prepare linear algebra variables for parallel solvers
 */
template<class LinearSolverTraits, class ParallelTraits,
         class Matrix, class Vector, class ParallelHelper>
void prepareLinearAlgebraParallel(Matrix& A, Vector& b, ParallelHelper& pHelper)
{
    if constexpr (ParallelTraits::isNonOverlapping)
    {
        // extend the matrix pattern such that it is usable for a parallel solver
        // and make right-hand side consistent
        using GridView = typename LinearSolverTraits::GridView;
        using DofMapper = typename LinearSolverTraits::DofMapper;
        static constexpr int dofCodim = LinearSolverTraits::dofCodim;
        ParallelMatrixHelper<Matrix, GridView, DofMapper, dofCodim> matrixHelper(pHelper.gridView(), pHelper.dofMapper());
        matrixHelper.extendMatrix(A, [&pHelper](auto idx){ return pHelper.isGhost(idx); });
        matrixHelper.sumEntries(A);

        ParallelVectorHelper<GridView, DofMapper, dofCodim> vectorHelper(pHelper.gridView(), pHelper.dofMapper());
        vectorHelper.makeNonOverlappingConsistent(b);
    }
}

/*!
 * \brief Prepare linear algebra variables for parallel solvers
 */
template<class LinearSolverTraits, class ParallelTraits,
         class Matrix, class Vector, class ParallelHelper>
void prepareLinearAlgebraParallel(Matrix& A, Vector& b,
                                  std::shared_ptr<typename ParallelTraits::Comm>& comm,
                                  std::shared_ptr<typename ParallelTraits::LinearOperator>& fop,
                                  std::shared_ptr<typename ParallelTraits::ScalarProduct>& sp,
                                  ParallelHelper& pHelper)
{
    prepareLinearAlgebraParallel<LinearSolverTraits, ParallelTraits>(A, b, pHelper);
    const auto category = ParallelTraits::isNonOverlapping ?
        Dune::SolverCategory::nonoverlapping : Dune::SolverCategory::overlapping;

    comm = std::make_shared<typename ParallelTraits::Comm>(pHelper.gridView().comm(), category);
    pHelper.createParallelIndexSet(*comm);
    fop = std::make_shared<typename ParallelTraits::LinearOperator>(A, *comm);
    sp = std::make_shared<typename ParallelTraits::ScalarProduct>(*comm);
}

} // end namespace Dumux

#endif // HAVE_MPI
#endif
