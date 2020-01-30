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
#ifndef DUMUX_PARALLELHELPERS_HH
#define DUMUX_PARALLELHELPERS_HH

#include <dune/common/version.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A parallel helper class providing a nonoverlapping
 *        decomposition of all degrees of freedom
 */
// operator that resets result to zero at constrained DOFS
template<class GridView, class LinearSolverTraits>
class ParallelISTLHelper
{
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

        bool fixedsize(int dim, int codim) const
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
     * \brief GatherScatter implementation that makes a right hand side in the box model consistent.
     */
    template<class V>
    class ConsistencyBoxGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<ConsistencyBoxGatherScatter<V>,typename V::block_type>
    {
    public:
        using DataType = typename V::block_type;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedsize;
        using BaseGatherScatter::size;

        ConsistencyBoxGatherScatter(V& container, const DofMapper& mapper)
        : BaseGatherScatter(mapper), container_(container)
        {}

        template<class MessageBuffer, class EntityType>
        void gather(MessageBuffer& buff, const EntityType& e) const
        {
            buff.write(container_[this->index(e)]);
        }

        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, std::size_t n)
        {
            typename V::block_type block;
            buff.read(block);
            container_[this->index(e)] += block;
        }
    private:
        V& container_;
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
        using BaseGatherScatter::fixedsize;
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
        using BaseGatherScatter::fixedsize;
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
        using BaseGatherScatter::fixedsize;
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
        using BaseGatherScatter::fixedsize;
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
        using BaseGatherScatter::fixedsize;
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
        : gridView_(gridView), mapper_(mapper), initialized_(false)
    {}

    [[deprecated("The verbose argument has no effect. Use ParallelISTLHelper(gridView, mapper) instead. Will be removed after 3.2!")]]
    ParallelISTLHelper(const GridView& gridView, const DofMapper& mapper, int verbose)
        : gridView_(gridView), mapper_(mapper), initialized_(false)
    {}

    // \brief Initializes the markers for ghosts and owners with the correct size and values.
    //
    void initGhostsAndOwners()
    {
        const auto rank = gridView_.comm().rank();
        isOwned_.resize(mapper_.size(), rank);
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

        initialized_ = true;
    }

    // keep only DOFs assigned to this processor
    template<class W>
    [[deprecated("This function has no effect. Will be removed after 3.2!")]]
    void mask(W& w) const
    {}

    // access to mask vector
    [[deprecated("Will be removed after 3.2!")]]
    std::size_t mask(std::size_t i) const
    { return isOwned_[i]; }

    // access to ghost vector
    [[deprecated("Will be removed after 3.2!")]]
    std::size_t ghost(std::size_t i) const
    { return isGhost_[i]; }

    bool isGhost(std::size_t i) const
    { return isGhost_[i] == ghostMarker_; }

    // \brief Make a vector of the box model consistent.
    template<class B, class A>
    void makeNonOverlappingConsistent(Dune::BlockVector<B,A>& v)
    {
#if HAVE_MPI
        ConsistencyBoxGatherScatter<Dune::BlockVector<B,A> > gs(v, mapper_);
        if (gridView_.comm().size() > 1)
            gridView_.communicate(gs, Dune::InteriorBorder_InteriorBorder_Interface,
                                  Dune::ForwardCommunication);
#endif
    }


#if HAVE_MPI

    template<typename MatrixType, typename Comm>
    [[deprecated("Use createParallelIndexSet(comm) instead. Will be removed after 3.2!")]]
    void createIndexSetAndProjectForAMG(MatrixType& m, Comm& c)
    { createParallelIndexSet(c); }

    /*!
     * \brief Creates a parallel index set
     *
     * \tparam Comm The type of the OwnerOverlapCopyCommunication
     * communicators.
     */
    template<class Comm>
    void createParallelIndexSet(Comm& comm)
    {
        if (!initialized_)
        {
            // This is the first time this function is called.
            // Therefore we need to initialize the marker vectors for ghosts and
            // owned dofs
            initGhostsAndOwners();
        }

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

#endif

    //! Return the dofMapper
    const DofMapper& dofMapper() const
    { return mapper_; }

    //! Return the gridView
    const GridView& gridView() const
    { return gridView_; }

    template<class Comm>
    Dune::OwnerOverlapCopyAttributeSet::AttributeSet
    getAttribute_(const Comm& comm, const bool isOwned, const bool isGhost) const
    {
        if (isOwned)
            return Dune::OwnerOverlapCopyAttributeSet::owner;
#if DUNE_VERSION_GTE(DUNE_ISTL, 2, 7)
        else if (isGhost && ( comm.category() ==
                                       static_cast<int>(Dune::SolverCategory::nonoverlapping)) )
#else
        else if (isGhost && ( comm.getSolverCategory() ==
                                       static_cast<int>(Dune::SolverCategory::nonoverlapping)) )
#endif
            return Dune::OwnerOverlapCopyAttributeSet::overlap;
        else
            return Dune::OwnerOverlapCopyAttributeSet::copy;
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

private:
    const GridView gridView_; //!< the grid view
    const DofMapper& mapper_; //!< the dof mapper
    std::vector<std::size_t> isOwned_; //!< vector to identify unique decomposition
    std::vector<std::size_t> isGhost_; //!< vector to identify ghost dofs
    bool initialized_; //!< whether isGhost and owner arrays are initialized

}; // class ParallelISTLHelper

/*!
 * \ingroup Linear
 * \brief Helper class for adding up matrix entries on border.
 * \tparam GridView The grid view to work on
 * \tparam LinearSolverTraits traits class
 */
template<class GridView, class LinearSolverTraits>
class EntityExchanger
{
    using Matrix = typename LinearSolverTraits::MType;
    enum { dim = GridView::dimension };
    enum { dofCodim = LinearSolverTraits::dofCodim };
    using Grid = typename GridView::Traits::Grid;
    using BlockType = typename Matrix::block_type;
    using IDS = typename Grid::Traits::GlobalIdSet;
    using IdType = typename IDS::IdType;
    using DofMapper = typename LinearSolverTraits::DofMapper;

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
     *  The MatrixPatternExchange class will find these entries and returns a vector "sparsity",
     *  that contains all missing connections.
     */
    struct MatrixPatternExchange
    : public Dune::CommDataHandleIF<MatrixPatternExchange, IdType>
    {
        //! Export type of data for message buffer
        using DataType = IdType;

        MatrixPatternExchange(const DofMapper& mapper,
                              const std::map<IdType,int>& globalToLocal,
                              const std::map<int,IdType>& localToGlobal, Matrix& A,
                              const ParallelISTLHelper<GridView, LinearSolverTraits>& helper)
        : mapper_(mapper), idToIndex_(globalToLocal), indexToID_(localToGlobal)
        , sparsity_(A.N()), A_(A), helper_(helper)
        {}

        /*!
         * \brief Returns true if data for given valid codim should be communicated
         */
        bool contains (int dim, int codim) const
        { return (codim == dofCodim); }

        /*!
         * \brief Returns true if size of data per entity of given dim and codim is a constant
         */
        bool fixedsize (int dim, int codim) const
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
                    if (!sparsity_[rowIdx].count(colIdx) && !helper_.isGhost(colIdx))
                        sparsity_[rowIdx].insert(colIdx);
                }
            }
        }

        /*!
         * \brief Get the communicated sparsity pattern
         * \return the vector with the sparsity pattern
         */
        std::vector<std::set<int>>& sparsity()
        { return sparsity_; }

    private:
        const DofMapper& mapper_;
        const std::map<IdType,int>& idToIndex_;
        const std::map<int,IdType>& indexToID_;
        std::vector<std::set<int> > sparsity_;
        Matrix& A_;
        const ParallelISTLHelper<GridView, LinearSolverTraits>& helper_;

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
    : public Dune::CommDataHandleIF<MatrixEntryExchange,MatrixEntry>
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

        /*!
         * \brief Returns true if size of data per entity of given dim and codim is a constant
         */
        bool fixedsize(int dim, int codim) const
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

    EntityExchanger(const GridView& gridView, const DofMapper& mapper)
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

    [[deprecated("Use extendMatrix instead. Will be removed after 3.2!")]]
    void getExtendedMatrix (Matrix& A, const ParallelISTLHelper<GridView, LinearSolverTraits>& helper)
    { extendMatrix(A, helper); }

    /*!
     * \brief communicates values for the sparsity pattern of the new matrix.
     * \param A Matrix to operate on.
     * \param helper ParallelelISTLHelper.
     */
    void extendMatrix(Matrix& A, const ParallelISTLHelper<GridView, LinearSolverTraits>& helper)
    {
        if (gridView_.comm().size() <= 1)
            return;

        Matrix tmp(A);
        std::size_t nnz = 0;
        // get entries from other processes
        MatrixPatternExchange datahandle(mapper_, idToIndex_, indexToID_, A, helper);
        gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface,
                                          Dune::ForwardCommunication);
        std::vector<std::set<int>>& sparsity = datahandle.sparsity();
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

}; // class EntityExchanger

#if HAVE_MPI
/*!
 * \brief Prepare linear algebra variables for parallel solvers
 */
template<class LinearSolverTraits, class Matrix, class Vector, class ParallelHelper>
void prepareLinearAlgebraParallel(Matrix& A, Vector& b,
                                  std::shared_ptr<typename LinearSolverTraits::Comm>& comm,
                                  std::shared_ptr<typename LinearSolverTraits::LinearOperator>& fop,
                                  std::shared_ptr<typename LinearSolverTraits::ScalarProduct>& sp,
                                  ParallelHelper& pHelper,
                                  const bool firstCall)
{
    const auto category = LinearSolverTraits::isNonOverlapping ?
                          Dune::SolverCategory::nonoverlapping
                        : Dune::SolverCategory::overlapping;

    if (LinearSolverTraits::isNonOverlapping && firstCall)
        pHelper.initGhostsAndOwners();

    comm = std::make_shared<typename LinearSolverTraits::Comm>(pHelper.gridView().comm(), category);

    if (LinearSolverTraits::isNonOverlapping)
    {
        // extend the matrix pattern such that it is usable for a parallel solver
        using GridView = std::decay_t<decltype(pHelper.gridView())>;
        EntityExchanger<GridView, LinearSolverTraits> exchanger(pHelper.gridView(), pHelper.dofMapper());
        exchanger.extendMatrix(A, pHelper);
        exchanger.sumEntries(A);
    }
    pHelper.createParallelIndexSet(*comm);

    fop = std::make_shared<typename LinearSolverTraits::LinearOperator>(A, *comm);
    sp = std::make_shared<typename LinearSolverTraits::ScalarProduct>(*comm);

    // make rhs consistent
    if (LinearSolverTraits::isNonOverlapping)
        pHelper.makeNonOverlappingConsistent(b);
}
#endif // HAVE_MPI

/*!
 * \brief Prepare linear algebra variables for sequential solvers
 */
template<class LinearSolverTraits, class Matrix>
void prepareLinearAlgebraSequential(Matrix& A,
                                    std::shared_ptr<typename LinearSolverTraits::Comm>& comm,
                                    std::shared_ptr<typename LinearSolverTraits::LinearOperator>& fop,
                                    std::shared_ptr<typename LinearSolverTraits::ScalarProduct>& sp)
{
    comm = std::make_shared<typename LinearSolverTraits::Comm>();
    fop = std::make_shared<typename LinearSolverTraits::LinearOperator>(A);
    sp = std::make_shared<typename LinearSolverTraits::ScalarProduct>();
}


/*!
 * \brief Prepare the linear algebra member variables.
 *
 * At compile time, correct constructor calls have to be chosen,
 * depending on whether the setting is parallel or sequential.
 * Since several template parameters are present, this cannot be solved
 * by a full function template specialization. Instead, class template
 * specialization has to be used.
 * This adapts example 4 from http://www.gotw.ca/publications/mill17.htm.
 *
 * This class template implements the function for the sequential case.
 *
 * \tparam isParallel decides if the setting is parallel or sequential
 */
template<class GridView, class AmgTraits, bool isParallel>
struct LinearAlgebraPreparator // TODO: Deprecating the whole struct always triggers a warning even it is not used. Remove after 3.2
{
    using ParallelHelper = ParallelISTLHelper<GridView, AmgTraits>;
    using Comm = typename AmgTraits::Comm;
    using LinearOperator = typename AmgTraits::LinearOperator;
    using ScalarProduct = typename AmgTraits::ScalarProduct;

    template<class Matrix, class Vector>
    [[deprecated("Use free function prepareLinearAlgebraSequential instead. Will be removed after 3.2!")]]
    static void prepareLinearAlgebra(Matrix& A, Vector& b,
                                     int& rank,
                                     std::shared_ptr<Comm>& comm,
                                     std::shared_ptr<LinearOperator>& fop,
                                     std::shared_ptr<ScalarProduct>& sp,
                                     ParallelHelper& pHelper,
                                     const bool firstCall)
    {
        prepareLinearAlgebraSequential<AmgTraits>(A, comm, fop, sp);
    }
};

#if HAVE_MPI
/*!
 * \brief Specialization for the parallel case.
 */
template<class GridView, class AmgTraits>
struct LinearAlgebraPreparator<GridView, AmgTraits, true> // TODO: Deprecating the whole struct always triggers a warning even it is not used. Remove after 3.2
{
    using ParallelHelper = ParallelISTLHelper<GridView, AmgTraits>;
    using Comm = typename AmgTraits::Comm;
    using LinearOperator = typename AmgTraits::LinearOperator;
    using ScalarProduct = typename AmgTraits::ScalarProduct;

    template<class Matrix, class Vector>
    [[deprecated("Use free function prepareLinearAlgebraParallel instead. Will be removed after 3.2!")]]
    static void prepareLinearAlgebra(Matrix& A, Vector& b,
                                     int& rank,
                                     std::shared_ptr<Comm>& comm,
                                     std::shared_ptr<LinearOperator>& fop,
                                     std::shared_ptr<ScalarProduct>& sp,
                                     ParallelHelper& pHelper,
                                     const bool firstCall)
    {
        prepareLinearAlgebraParallel<AmgTraits>(A, b, rank, comm, fop, sp, pHelper, firstCall);
        rank = comm->communicator().rank();
    }

}; // parallel LinearAlgebraPreparator

#endif // HAVE_MPI

} // end namespace Dumux

#endif // DUMUX_PARALLELHELPERS_HH
