// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Provides a helper class for nonoverlapping decomposition
 */
#ifndef DUMUX_LINEAR_PARALLELHELPERS_HH
#define DUMUX_LINEAR_PARALLELHELPERS_HH

#include <dune/common/exceptions.hh>
#include <dune/geometry/dimension.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>
#include <dumux/common/gridcapabilities.hh>

namespace Dumux::Detail {

template<class LinearSolverTraits, bool canCommunicate = false>
class ParallelISTLHelperImpl {};

template<class LinearSolverTraits>
class ParallelISTLHelperImpl<LinearSolverTraits, true>
{
    using GridView = typename LinearSolverTraits::GridView;
    using DofMapper = typename LinearSolverTraits::DofMapper;
    static constexpr int dofCodim = LinearSolverTraits::dofCodim;

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
    : public BaseGatherScatter
    , public Dune::CommDataHandleIF<GhostGatherScatter, std::size_t>
    {
    public:
        using DataType = std::size_t;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedSize;
        using BaseGatherScatter::size;

        GhostGatherScatter(std::vector<std::size_t>& ranks, const DofMapper& mapper)
        : BaseGatherScatter(mapper)
        , ranks_(ranks)
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
    : public BaseGatherScatter
    , public Dune::CommDataHandleIF<InteriorBorderGatherScatter, std::size_t>
    {
    public:
        using DataType = std::size_t;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedSize;
        using BaseGatherScatter::size;

        InteriorBorderGatherScatter(std::vector<std::size_t>& ranks, const DofMapper& mapper)
        : BaseGatherScatter(mapper)
        , ranks_(ranks)
        {}

        /*!
         * \brief Pack data (to be send to other processes) from this process to message buffer
         */
        template<class MessageBuffer, class EntityType>
        void gather(MessageBuffer& buff, const EntityType& e) const
        {
            auto& data = ranks_[this->index(e)];
            if (this->isNeitherInteriorNorBorderEntity(e))
                data = ghostMarker_;
            buff.write(data);
        }

        /*!
         * \brief Unpack data (received from other process) from message buffer to this process
         */
        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, std::size_t n)
        {
            std::size_t x;
            buff.read(x);
            auto& data = ranks_[this->index(e)];

            // we leave ghost unchanged
            // for other dofs, the process with the lowest rank
            // is assigned to be the (unique) owner
            using std::min;
            data = this->isNeitherInteriorNorBorderEntity(e) ? x : min(data, x);
        }
    private:
        std::vector<std::size_t>& ranks_;
    };

    /*!
     * \brief GatherScatter handle for finding out about neighbouring processor ranks.
     *
     */
    struct NeighbourGatherScatter
    : public BaseGatherScatter
    , public Dune::CommDataHandleIF<NeighbourGatherScatter, int>
    {
        using DataType = int;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedSize;
        using BaseGatherScatter::size;

        NeighbourGatherScatter(const DofMapper& mapper, int rank, std::set<int>& neighbours)
        : BaseGatherScatter(mapper)
        , rank_(rank)
        , neighbours_(neighbours)
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
    : public BaseGatherScatter
    , public Dune::CommDataHandleIF<SharedGatherScatter, int>
    {
        using DataType = int;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedSize;
        using BaseGatherScatter::size;

        SharedGatherScatter(std::vector<int>& shared, const DofMapper& mapper)
        : BaseGatherScatter(mapper)
        , shared_(shared)
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
    : public BaseGatherScatter
    , public Dune::CommDataHandleIF<GlobalIndexGatherScatter<GlobalIndex>, GlobalIndex>
    {
        using DataType = GlobalIndex;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedSize;
        using BaseGatherScatter::size;

        GlobalIndexGatherScatter(std::vector<GlobalIndex>& globalIndices, const DofMapper& mapper)
        : BaseGatherScatter(mapper)
        , globalIndices_(globalIndices)
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
            DataType& data = globalIndices_[this->index(e)];
            data = min(data, x);
        }
    private:
        std::vector<GlobalIndex>& globalIndices_;
    };

public:

    ParallelISTLHelperImpl(const GridView& gridView, const DofMapper& mapper)
    : gridView_(gridView), mapper_(mapper)
    {
        if constexpr (Detail::canCommunicate<typename GridView::Traits::Grid, dofCodim>)
            initGhostsAndOwners_();
        else
            DUNE_THROW(Dune::InvalidStateException,
                "Cannot initialize parallel helper for a grid that cannot communicate codim-" << dofCodim << "-entities."
            );
    }

    bool isGhost(std::size_t i) const
    { return isGhost_[i] == ghostMarker_; }

    bool isOwned(std::size_t i) const
    { return isOwned_[i] == 1; }

    /*!
     * \brief Creates a parallel index set for a communicator
     * \param comm The OwnerOverlapCopyCommunication communicators
     *
     * The parallel index set contains for each dof index,
     * the triplet (globalIndex, (localIndex, attribute))
     * where attribute is owner, overlap or copy. Each dof
     * is uniquely owned by exactly one process. This allows
     * to apply local operators additively to get a global operator
     * without communication.
     */
    template<class Comm>
    void createParallelIndexSet(Comm& comm) const
    {
        if constexpr (Detail::canCommunicate<typename GridView::Traits::Grid, dofCodim>)
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

            // starting from start (offset), number global indices owned by this process consecutively
            std::vector<GlobalIndex> globalIndices(mapper_.size(), std::numeric_limits<GlobalIndex>::max());
            for (std::size_t i = 0; i < globalIndices.size(); ++i)
            {
                if (isOwned_[i] == 1 && isShared[i])
                {
                    globalIndices[i] = start;
                    ++start;
                }
            }

            // publish global indices for the shared DOFS to other processors.
            GlobalIndexGatherScatter<GlobalIndex> gigs(globalIndices, mapper_);
            gridView_.communicate(gigs, Dune::All_All_Interface, Dune::ForwardCommunication);

            // put the information into the parallel index set
            resizeAndFillIndexSet_(comm, globalIndices);

            // Compute neighbours using communication
            std::set<int> neighbours;
            NeighbourGatherScatter ngs(mapper_, gridView_.comm().rank(), neighbours);
            gridView_.communicate(ngs, Dune::All_All_Interface, Dune::ForwardCommunication);
            comm.remoteIndices().setNeighbours(neighbours);

            comm.remoteIndices().template rebuild<false>();
        }
        else
            DUNE_THROW(Dune::InvalidStateException,
                "Cannot build parallel index set for a grid that cannot communicate codim-" << dofCodim << "-entities."
            );
    }

    //! Return the dofMapper
    const DofMapper& dofMapper() const
    { return mapper_; }

    //! Return the gridView
    const GridView& gridView() const
    { return gridView_; }

private:
     /*!
     * \brief Initialize ghost and owner flags
     *
     * All dofs on ghost entities are marked.
     * All dofs are assigned a unique owner process. This will help us to create a unique partition
     * and unique representations of a vector according to Definition 2.5 in Blatt and Bastian (2009)
     * https://doi.org/10.1504/IJCSE.2008.021112
     */
    void initGhostsAndOwners_()
    {
        const auto rank = gridView_.comm().rank();
        isOwned_.assign(mapper_.size(), rank);

        // find out about ghosts
        GhostGatherScatter ggs(isOwned_, mapper_);

        if (gridView_.comm().size() > 1)
            gridView_.communicate(ggs, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);

        // isGhost_ contains rank when the dof is not on a ghost entity (interior or border)
        // and ghostMarker_ when the dof is on a ghost
        isGhost_ = isOwned_;

        // partition interior/border uniquely
        // after the communication each dof is assigned to one unique owner process rank
        InteriorBorderGatherScatter dh(isOwned_, mapper_);

        if (gridView_.comm().size() > 1)
            gridView_.communicate(dh, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);

        // convert vector into mask vector
        for (auto& v : isOwned_)
            v = (v == rank) ? 1 : 0;
    }

    template<class Comm, class GlobalIndices>
    void resizeAndFillIndexSet_(Comm& comm, const GlobalIndices& globalIndices) const
    {
        comm.indexSet().beginResize();

        // add triplets characterizing each dof in the parallel index set
        // (globalIndex, (localIndex, attribute)) where attribute is owner, overlap or copy
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
        if (isOwned) // this process owns the dof (uniquely)
            return Dune::OwnerOverlapCopyAttributeSet::owner;
        else if (isGhost && (comm.category() == Dune::SolverCategory::nonoverlapping) )
            return Dune::OwnerOverlapCopyAttributeSet::overlap;
        else
            return Dune::OwnerOverlapCopyAttributeSet::copy;
    }

    const GridView gridView_; //!< the grid view
    const DofMapper& mapper_; //!< the dof mapper
    //! vector to identify unique decomposition (1: this rank owns the dof; 0: owned by remote process)
    std::vector<std::size_t> isOwned_;
    //!< vector to identify ghost dofs (ghostMarker_: this dof is on a ghost entity; contains and unspecified value otherwise)
    std::vector<std::size_t> isGhost_;

};

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief A parallel helper class providing a parallel
 *        decomposition of all degrees of freedom
 */
template<class LinearSolverTraits>
using ParallelISTLHelper =
    Detail::ParallelISTLHelperImpl<
        LinearSolverTraits, LinearSolverTraits::canCommunicate
    >;


template<class GridView, class DofMapper, int dofCodim>
class ParallelVectorHelper
{
public:
    ParallelVectorHelper(const GridView& gridView, const DofMapper& mapper)
    : gridView_(gridView), mapper_(mapper)
    {}

    //! \brief Make a vector consistent for non-overlapping domain decomposition methods
    template<class Block, class Alloc>
    void makeNonOverlappingConsistent(Dune::BlockVector<Block, Alloc>& v) const
    {
        if constexpr (Detail::canCommunicate<typename GridView::Traits::Grid, dofCodim>)
        {
            VectorCommDataHandleSum<DofMapper, Dune::BlockVector<Block, Alloc>, dofCodim, Block> gs(mapper_, v);
            if (gridView_.comm().size() > 1)
                gridView_.communicate(gs, Dune::InteriorBorder_InteriorBorder_Interface,
                                    Dune::ForwardCommunication);
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Cannot call makeNonOverlappingConsistent for a grid that cannot communicate codim-" << dofCodim << "-entities.");
    }

    //! \brief Make a vector consistent for non-overlapping domain decomposition methods
    template<class... Blocks>
    void makeNonOverlappingConsistent(Dune::MultiTypeBlockVector<Blocks...>& v) const
    {
        DUNE_THROW(Dune::NotImplemented, "makeNonOverlappingConsistent for Dune::MultiTypeBlockVector");
    }

private:
    const GridView gridView_; //!< the grid view
    const DofMapper& mapper_; //!< the dof mapper
};

/*!
 * \ingroup Linear
 * \brief Helper class for adding up matrix entries for border entities
 *
 * Border means all degrees of freedom located
 * on lower-dimensional entities (faces, edges, vertices)
 * that form the the process boundary
 */
template<class Matrix, class GridView,
         class RowDofMapper, int rowDofCodim>
class ParallelMatrixHelper
{
    using IdType = typename GridView::Traits::Grid::Traits::GlobalIdSet::IdType;

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
    template<class ColIsGhostFunc>
    struct MatrixPatternExchange
    : public Dune::CommDataHandleIF<MatrixPatternExchange<ColIsGhostFunc>, IdType>
    {
        //! Export type of data for message buffer
        using DataType = IdType;

        /*!
         * \brief Construct a new sparsity pattern exchange handle
         *
         * \param rowEntityMapper an index mapper for entities associated with row entities
         * \param globalToLocal map (column dof) from global ID to processor-local index
         * \param localToGlobal map (column dof) from  processor-local index to global ID
         * \param A the matrix for which we want to extend the pattern
         * \param[out] sparsityPattern the extended sparsity pattern
         * \param isGhostColumDof a function that returns whether a given column index is associated with a ghost entity
         */
        MatrixPatternExchange(const RowDofMapper& rowEntityMapper,
                              const std::map<IdType,int>& globalToLocal,
                              const std::map<int,IdType>& localToGlobal,
                              Matrix& A,
                              std::vector<std::set<int>>& sparsityPattern,
                              const ColIsGhostFunc& isGhostColumDof)
        : rowEntityMapper_(rowEntityMapper), idToIndex_(globalToLocal), indexToID_(localToGlobal)
        , sparsityPattern_(sparsityPattern), A_(A), isGhostColumDof_(isGhostColumDof)
        {
            sparsityPattern_.resize(A.N());
        }

        /*!
         * \brief Returns true if data for given valid codim should be communicated
         */
        bool contains (int dim, int codim) const
        { return (codim == rowDofCodim); }

        //! returns true if size per entity of given dim and codim is a constant
        bool fixedSize(int dim, int codim) const
        { return false; }

        /*!
         * \brief How many objects of type DataType have to be sent for a given entity
         */
        template<class EntityType>
        std::size_t size(EntityType& e) const
        {
            const auto rowIdx = rowEntityMapper_.index(e);

            // all column entity indices of this row that are in the index set
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
            const auto rowIdx = rowEntityMapper_.index(e);

            // send all column entity ids of this row that are in the index set
            for (auto colIt = A_[rowIdx].begin(); colIt != A_[rowIdx].end(); ++colIt)
                if (auto it = indexToID_.find(colIt.index()); it != indexToID_.end())
                    buff.write(it->second);
        }

        /*!
         * \brief Unpack data from message buffer to user
         */
        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, std::size_t n)
        {
            const auto rowIdx = rowEntityMapper_.index(e);
            for (std::size_t k = 0; k < n; k++)
            {
                // receive n column entity IDs
                IdType id;
                buff.read(id);

                // only add entries that are contained in the index set
                if (const auto it = idToIndex_.find(id); it != idToIndex_.end())
                {
                    const auto colIdx = it->second; // get local column index
                    // add this entry (if it doesn't exist yet)
                    // (and only if the column dof is not associated with a ghost entity)
                    if (!isGhostColumDof_(colIdx))
                        // std::set takes care that each index only is inserted once
                        sparsityPattern_[rowIdx].insert(colIdx);
                }
            }
        }

    private:
        const RowDofMapper& rowEntityMapper_;
        const std::map<IdType, int>& idToIndex_;
        const std::map<int, IdType>& indexToID_;
        std::vector<std::set<int>>& sparsityPattern_;
        Matrix& A_;
        const ColIsGhostFunc& isGhostColumDof_;

    }; // class MatrixPatternExchange

    //! Local matrix blocks associated with the global id set
    struct MatrixEntry
    {
        IdType first;
        typename Matrix::block_type second;
    };

    //! A DataHandle class to exchange matrix entries
    struct MatrixEntryExchange
    : public Dune::CommDataHandleIF<MatrixEntryExchange, MatrixEntry>
    {
        //! Export type of data for message buffer
        using DataType = MatrixEntry;

        MatrixEntryExchange(const RowDofMapper& rowEntityMapper,
                            const std::map<IdType, int>& globalToLocal,
                            const std::map<int, IdType>& localToGlobal,
                            Matrix& A)
        : rowEntityMapper_(rowEntityMapper), idToIndex_(globalToLocal), indexToID_(localToGlobal), A_(A)
        {}

        /*!
         * \brief Returns true if data for given valid codim should be communicated
         */
        bool contains(int dim, int codim) const
        { return (codim == rowDofCodim); }

        //! returns true if size per entity of given dim and codim is a constant
        bool fixedSize(int dim, int codim) const
        { return false; }

        /*!
         * \brief How many objects of type DataType have to be sent for a given entity
         */
        template<class EntityType>
        std::size_t size(EntityType& e) const
        {
            const auto rowIdx = rowEntityMapper_.index(e);
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
            const auto rowIdx = rowEntityMapper_.index(e);
            // send all matrix entries for which the column index is in the index set
            for (auto colIt = A_[rowIdx].begin(); colIt != A_[rowIdx].end(); ++colIt)
                if (auto it = indexToID_.find(colIt.index()); it != indexToID_.end())
                    buff.write(MatrixEntry{ it->second, *colIt });
        }

        /*!
         * \brief Unpack data from message buffer to user
         */
        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, std::size_t n)
        {
            const auto rowIdx = rowEntityMapper_.index(e);
            for (std::size_t k = 0; k < n; k++)
            {
                MatrixEntry m;
                buff.read(m);
                const auto& [colDofID, matrixBlock] = m; // unpack

                // only add entries in the index set and for which A has allocated memory
                if (auto it = idToIndex_.find(colDofID); it != idToIndex_.end())
                    if (A_[rowIdx].find(it->second) != A_[rowIdx].end())
                        A_[rowIdx][it->second] += matrixBlock;
            }
        }

    private:
        const RowDofMapper& rowEntityMapper_;
        const std::map<IdType, int>& idToIndex_;
        const std::map<int, IdType>& indexToID_;
        Matrix& A_;

    }; // class MatrixEntryExchange

public:

    ParallelMatrixHelper(const GridView& gridView, const RowDofMapper& mapper)
    : gridView_(gridView), mapper_(mapper)
    {
        idToIndex_.clear();
        indexToID_.clear();

        std::vector<bool> handledDof(gridView_.size(rowDofCodim), false);
        for (const auto& element : elements(gridView_))
        {
            for (int i = 0; i < element.subEntities(rowDofCodim); ++i)
            {
                const auto entity = element.template subEntity<rowDofCodim>(i);
                if (entity.partitionType() == Dune::BorderEntity)
                {
                    const auto localRowIdx = mapper_.index(entity);
                    if (!handledDof[localRowIdx])
                    {
                        IdType dofIdxGlobal = gridView_.grid().globalIdSet().id(entity);
                        idToIndex_.emplace(dofIdxGlobal, localRowIdx);
                        indexToID_.emplace(localRowIdx, dofIdxGlobal);
                    }
                }
            }
        }
    }

    /*!
     * \brief communicates values for the sparsity pattern of the new matrix.
     * \param A Matrix to operate on.
     * \param isGhost Function returning if a column dof index is on a ghost entity
     */
    template<class IsGhostFunc>
    void extendMatrix(Matrix& A, const IsGhostFunc& isGhost)
    {
        if constexpr (Detail::canCommunicate<typename GridView::Grid, rowDofCodim>)
        {
            if (gridView_.comm().size() <= 1)
                return;

            // make a copy of the matrix as originally assembled
            Matrix matrixAsAssembled(A);

            // first get entries in sparsity pattern from other processes
            std::size_t numNonZeroEntries = 0;
            std::vector<std::set<int>> sparsityPattern; // column indices for every row
            MatrixPatternExchange<IsGhostFunc> dataHandle(mapper_, idToIndex_, indexToID_, A, sparsityPattern, isGhost);
            gridView_.communicate(
                dataHandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication
            );

            // add own entries to the sparsity pattern and count number of non-zeros
            for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
            {
                const auto colEndIt = A[rowIt.index()].end();
                for (auto colIt = rowIt->begin(); colIt != colEndIt; ++colIt)
                    if (!sparsityPattern[rowIt.index()].count(colIt.index()))
                        sparsityPattern[rowIt.index()].insert(colIt.index());

                numNonZeroEntries += sparsityPattern[rowIt.index()].size();
            }

            // insert any additional entries into the matrix
            A.setSize(matrixAsAssembled.N(), matrixAsAssembled.M(), numNonZeroEntries);
            A.setBuildMode(Matrix::row_wise);
            auto citer = A.createbegin();
            for (auto i = sparsityPattern.begin(), end = sparsityPattern.end(); i!=end; ++i, ++citer)
                for (auto si = i->begin(), send = i->end(); si!=send; ++si)
                    citer.insert(*si);

            // reset matrix to contain original values
            A = 0;
            const auto rowEndIt = matrixAsAssembled.end();
            for (auto rowIt = matrixAsAssembled.begin(); rowIt != rowEndIt; ++rowIt)
                for (auto colIt = matrixAsAssembled[rowIt.index()].begin(); colIt != matrixAsAssembled[rowIt.index()].end(); ++colIt)
                    A[rowIt.index()][colIt.index()] = *colIt;

            // The matrix A has now a possible extended sparsity pattern but any additional entries are
            // initialized to zero and have to be filled by the sumEntries function
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "Cannot call extendMatrix for a grid that cannot communicate codim-" << rowDofCodim << "-entities.");
    }

    /*!
     * \brief Sums up the entries corresponding to border entities (usually vertices or faces)
     * \param A Matrix to operate on
     *
     * The idea is as follows: (Blatt and Bastian (2009) https://doi.org/10.1504/IJCSE.2008.021112)
     * The local matrix operator stores for each row that corresponds to a dof (uniquely) owned by this process
     * the full row of the global operator with the same entries as the global operator.

     * This, together with some masking procedure (end of Section 2.4.1) allows to compute a matrix-vector product
     * given a consistent vector representation such that the result is in a valid representation. Each valid
     * representation can be transformed to a additive unique representation by setting all entries for
     * dofs that are not (uniquely) owned by the process to zero.
     */
    void sumEntries(Matrix& A)
    {
        if constexpr (Detail::canCommunicate<typename GridView::Grid, rowDofCodim>)
        {
            if (gridView_.comm().size() <= 1)
                return;

            MatrixEntryExchange dataHandle(mapper_, idToIndex_, indexToID_, A);
            gridView_.communicate(
                dataHandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication
            );
        }
        else
            DUNE_THROW(Dune::InvalidStateException,
                "Cannot call sumEntries for a grid that cannot communicate codim-" << rowDofCodim << "-entities."
            );
    }

private:
    const GridView gridView_;
    const RowDofMapper& mapper_;
    std::map<IdType, int> idToIndex_;
    std::map<int, IdType> indexToID_;

};

/*!
 * \brief Prepare a matrix for parallel solvers
 */
template<class LinearSolverTraits, class ParallelTraits,
         class Matrix, class ParallelHelper>
void prepareMatrixParallel(Matrix& A, ParallelHelper& pHelper)
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
    }
}

/*!
 * \brief Prepare a vector for parallel solvers
 */
template<class LinearSolverTraits, class ParallelTraits,
         class Vector, class ParallelHelper>
void prepareVectorParallel(Vector& b, ParallelHelper& pHelper)
{
    if constexpr (ParallelTraits::isNonOverlapping)
    {
        // extend the matrix pattern such that it is usable for a parallel solver
        // and make right-hand side consistent
        using GridView = typename LinearSolverTraits::GridView;
        using DofMapper = typename LinearSolverTraits::DofMapper;
        static constexpr int dofCodim = LinearSolverTraits::dofCodim;
        ParallelVectorHelper<GridView, DofMapper, dofCodim> vectorHelper(pHelper.gridView(), pHelper.dofMapper());
        vectorHelper.makeNonOverlappingConsistent(b);
    }
}

/*!
 * \brief Prepare linear algebra variables for parallel solvers
 */
template<class LinearSolverTraits, class ParallelTraits,
         class Matrix, class Vector, class ParallelHelper>
void prepareLinearAlgebraParallel(Matrix& A, Vector& b, ParallelHelper& pHelper)
{
    prepareMatrixParallel<LinearSolverTraits, ParallelTraits>(A, pHelper);
    prepareVectorParallel<LinearSolverTraits, ParallelTraits>(b, pHelper);
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

#endif
