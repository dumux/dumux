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
        size_t size(EntityType& e) const
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
        void scatter(MessageBuffer& buff, const EntityType& e, size_t n)
        {
            typename V::block_type block;
            buff.read(block);
            container_[this->index(e)] += block;
        }
    private:
        V& container_;
    };


    /*!
     * \brief Writes 1<<24 to each data item (of the container) that is gathered or scattered
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
            std::size_t& data = ranks_[this->index(e)]; // TODO: reference to size_t ?
            if (this->isNeitherInteriorNorBorderEntity(e))
                data = (1<<24);
            buff.write(data);
        }

        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, size_t n)
        {
            std::size_t x;
            std::size_t& data = ranks_[this->index(e)];
            buff.read(x); // TODO what is this good for?
            if (this->isNeitherInteriorNorBorderEntity(e))
                data = (1<<24);
        }
    private:
        std::vector<std::size_t>& ranks_;
    };

    /*!
     * \brief GatherScatter handle that sets 1<<24 for data items neither associated to
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
            std::size_t& data = ranks_[this->index(e)]; // TODO: reference to size_t ?
            if (this->isNeitherInteriorNorBorderEntity(e))
                data = (1<<24);
            buff.write(data);
        }

        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, std::size_t n)
        {
            using std::min;
            std::size_t x;
            std::size_t& data = ranks_[this->index(e)];
            buff.read(x);
            if (this->isNeitherInteriorNorBorderEntity(e))
                data = x;
            else
                data = min(data,x);
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
        using DataType = int; // TODO: size_t?
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
        void scatter(MessageBuffer& buff, const EntityType& e, size_t n)
        {
            int x; // TODO: size_t?
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
        void scatter(MessageBuffer& buff, const EntityType &e, size_t n)
        {
            int x;
            buff.read(x);
            int& data = shared_[this->index(e)];
            data = data || x;
        }
    private:
        std::vector<int>& shared_;
    };

    /*!
     * \brief GatherScatter handle for finding out about neighbouring processor ranks.
     *
     */
    template<typename GI>
    struct GlobalIndexGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<GlobalIndexGatherScatter<GI>, GI>
    {
        using DataType = GI;
        using BaseGatherScatter::contains;
        using BaseGatherScatter::fixedsize;
        using BaseGatherScatter::size;

        GlobalIndexGatherScatter(std::vector<GI>& gindices, const DofMapper& mapper)
        : BaseGatherScatter(mapper), gindices_(gindices)
        {}

        template<class MessageBuffer, class EntityType>
        void gather(MessageBuffer& buff, const EntityType& e) const
        {
            buff.write(gindices_[this->index(e)]);
        }

        template<class MessageBuffer, class EntityType>
        void scatter(MessageBuffer& buff, const EntityType& e, size_t n)
        {
            using std::min;
            DataType x;
            buff.read(x);
            gindices_[this->index(e)] = min(gindices_[this->index(e)], x);
        }
    private:
        std::vector<GI>& gindices_;
    };

public:

    ParallelISTLHelper (const GridView& gridView, const DofMapper& mapper, int verbose = 1)
        : gridView_(gridView), mapper_(mapper), verbose_(verbose), initialized_(false)
    {}

    // \brief Initializes the markers for ghosts and owners with the correct size and values.
    //
    void initGhostsAndOwners()
    {
        owner_.resize(mapper_.size(),
                      gridView_.comm().rank());
        isGhost_.resize(mapper_.size(),0.0);
        // find out about ghosts
        GhostGatherScatter ggs(owner_, mapper_);

        if (gridView_.comm().size() > 1)
            gridView_.communicate(ggs, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);

        isGhost_ = owner_;

        // partition interior/border
        InteriorBorderGatherScatter dh(owner_, mapper_);

        if (gridView_.comm().size() > 1)
            gridView_.communicate(dh, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);

        // convert vector into mask vector
        for (auto& v : owner_)
            v = (v == gridView_.comm().rank()) ? 1.0 : 0.0;

        initialized_ = true;
    }

    // keep only DOFs assigned to this processor
    template<class W>
    void mask(W& w) const
    {
        auto v1 = w.begin();
        for (auto v2 = owner_.begin(), vend = owner_.end(); v2 != vend; ++v1, ++v2)
            v1 *= v2;
    }

    // access to mask vector
    double mask(std::size_t i) const
    { return owner_[i]; }

    // access to ghost vector
    std::size_t ghost (std::size_t i) const
    { return isGhost_[i]; }

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

    /*!
     * \brief Creates a matrix suitable for parallel AMG and the parallel information
     *
     * \tparam MatrixType The type of the ISTL matrix used.
     * \tparam Comm The type of the OwnerOverlapCopyCommunication
     * \param m The local matrix.
     * \param c The parallel information object providing index set, interfaces and
     * communicators.
     */
    template<typename MatrixType, typename Comm>
    [[deprecated("Use createParallelIndexSet(comm) instead. Will be removed after 3.1!")]]
    void createIndexSetAndProjectForAMG(MatrixType& m, Comm& c)
    { createParallelIndexSet(c); }

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

        // First find out which dofs we share with other processors
        std::vector<int> sharedDofs(mapper_.size(), false);

        SharedGatherScatter sgs(sharedDofs, mapper_);

        if (gridView_.comm().size() > 1)
            gridView_.communicate(sgs, Dune::All_All_Interface,
                                  Dune::ForwardCommunication);

        // Count shared dofs that we own
        using GlobalIndex = typename Comm::ParallelIndexSet::GlobalIndex;
        GlobalIndex count = 0;
        auto owned = owner_.begin();

        for (auto v = sharedDofs.begin(), vend=sharedDofs.end(); v != vend; ++v, ++owned)
            if (*v && *owned == 1)
                ++count;

        Dune::dverb << gridView_.comm().rank() << ": shared count is " << count.touint()
                    << std::endl;

        std::vector<GlobalIndex> counts(gridView_.comm().size());
        gridView_.comm().allgather(&count, 1, &(counts[0]));

        // Compute start index start_p = \sum_{i=0}^{i<p} counts_i
        GlobalIndex start = 0;
        for (int i = 0; i < gridView_.comm().rank(); ++i)
            start = start + counts[i];

        std::vector<GlobalIndex> scalarIndices(mapper_.size(),
                                               std::numeric_limits<GlobalIndex>::max());

        auto shared = sharedDofs.begin();
        auto index = scalarIndices.begin();

        for (auto i = owner_.begin(), iend = owner_.end(); i != iend; ++i, ++shared, ++index)
        {
            if (*i==1 && *shared)
            {
                *index = start;
                ++start;
            }
        }

        // publish global indices for the shared DOFS to other processors.
        GlobalIndexGatherScatter<GlobalIndex> gigs(scalarIndices, mapper_);
        if (gridView_.comm().size() > 1)
            gridView_.communicate(gigs, Dune::All_All_Interface,
                                  Dune::ForwardCommunication);


        // Setup the index set
        comm.indexSet().beginResize();
        index = scalarIndices.begin();
        auto ghost = isGhost_.begin();

        for (auto i = owner_.begin(), iend = owner_.end(); i != iend; ++i, ++ghost, ++index)
        {
            Dune::OwnerOverlapCopyAttributeSet::AttributeSet attr;
            if (*index != std::numeric_limits<GlobalIndex>::max())
            {
                // global index exist in index set
                if (*i > 0)
                {
                    // This dof is managed by us.
                    attr = Dune::OwnerOverlapCopyAttributeSet::owner;
                }
    #if DUNE_VERSION_GTE(DUNE_ISTL, 2, 7)
                else if ( *ghost==(1<<24) && ( comm.category() ==
                                               static_cast<int>(Dune::SolverCategory::nonoverlapping)) )
    #else
                else if ( *ghost==(1<<24) && ( comm.getSolverCategory() ==
                                               static_cast<int>(Dune::SolverCategory::nonoverlapping)) )
    #endif
                {
                    //use attribute overlap for ghosts in novlp grids
                    attr = Dune::OwnerOverlapCopyAttributeSet::overlap;
                }
                else
                {
                    attr = Dune::OwnerOverlapCopyAttributeSet::copy;
                }
                comm.indexSet().add(*index, typename Comm::ParallelIndexSet::LocalIndex(i-owner_.begin(), attr));
            }
        }
        comm.indexSet().endResize();

        // Compute neighbours using communication
        std::set<int> neighbours;
        NeighbourGatherScatter ngs(mapper_, gridView_.comm().rank(),
                                   neighbours);

        if (gridView_.comm().size() > 1)
            gridView_.communicate(ngs, Dune::All_All_Interface,
                                  Dune::ForwardCommunication);

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

private:
    const GridView gridView_; //!< the grid view
    const DofMapper& mapper_; //!< the dof mapper
    std::vector<std::size_t> owner_; //!< vector to identify unique decomposition
    std::vector<std::size_t> isGhost_; //!< vector to identify ghost dofs
    int verbose_; //!< verbosity
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
    using RowIterator = typename Matrix::RowIterator;
    using ColIterator = typename Matrix::ColIterator;
    using DofMapper = typename LinearSolverTraits::DofMapper;

public:
    /*! \brief Constructor. Sets up the local to global relations.
      \param[in] gridView The gridView on which we are operating
      \param[in] mapper The local dof mapper
    */
    EntityExchanger(const GridView& gridView, const DofMapper& mapper)
    : gridView_(gridView), mapper_(mapper)
    {
        gid2Index_.clear();
        index2GID_.clear();

        for (const auto& entity : entities(gridView_, Dune::Codim<dofCodim>()))
        {
            if (entity.partitionType() == Dune::BorderEntity)
            {
                int localIdx = mapper_.index(entity);
                IdType dofIdxGlobal = gridView_.grid().globalIdSet().id(entity);

                std::pair<IdType, int> g2iPair(dofIdxGlobal, localIdx);
                gid2Index_.insert(g2iPair);

                std::pair<int, IdType> i2gPair(localIdx, dofIdxGlobal);
                index2GID_.insert(i2gPair);
            }
        }
    }

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
     *  The MatPatternExchange class will find these entries and returns a vector "sparsity",
     *  that contains all missing connections.
     */
    class MatPatternExchange
        : public Dune::CommDataHandleIF<MatPatternExchange, IdType>
    {
        using RowIterator = typename Matrix::RowIterator;
        using ColIterator = typename Matrix::ColIterator;
    public:
        //! Export type of data for message buffer
        using DataType = IdType;

        /** \brief Constructor
            \param[in] mapper The local dof mapper.
            \param[in] g2i Global to local index map.
            \param[in] i2g Local to global index map.
            \param[in] A Matrix to operate on.
            \param[in] helper parallel istl helper.
        */
        MatPatternExchange (const DofMapper& mapper,
                            const std::map<IdType,int>& g2i,
                            const std::map<int,IdType>& i2g, Matrix& A,
                            const ParallelISTLHelper<GridView, LinearSolverTraits>& helper)
            : mapper_(mapper), gid2Index_(g2i), index2GID_(i2g),
              sparsity_(A.N()), A_(A), helper_(helper)
        {}

        /*!
         * \brief Returns true if data for given valid codim should be communicated
         */
        bool contains (int dim, int codim) const
        {
            return (codim==dofCodim);
        }

        /*!
         * \brief Returns true if size of data per entity of given dim and codim is a constant
         */
        bool fixedsize (int dim, int codim) const
        {
            return false;
        }

        /*!
         * \brief How many objects of type DataType have to be sent for a given entity
         */
        template<class EntityType>
        size_t size (EntityType& e) const
        {
            int i = mapper_.index(e);
            int n = 0;
            for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
                typename std::map<int,IdType>::const_iterator it = index2GID_.find(j.index());
                if (it != index2GID_.end())
                    n++;
            }

            return n;
        }

        /*!
         * \brief Pack data from user to message buffer
         */
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
            int i = mapper_.index(e);
            for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
                typename std::map<int,IdType>::const_iterator it=index2GID_.find(j.index());
                if (it != index2GID_.end())
                    buff.write(it->second);
            }

        }

        /*!
         * \brief Unpack data from message buffer to user
         */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
            int i = mapper_.index(e);
            for (size_t k = 0; k < n; k++)
            {
                IdType id;
                buff.read(id);
                // only add entries corresponding to border entities
                typename std::map<IdType,int>::const_iterator it = gid2Index_.find(id);
                if (it != gid2Index_.end()
                    && sparsity_[i].find(it->second) == sparsity_[i].end()
                    && helper_.ghost(it->second) != 1<<24)
                    sparsity_[i].insert(it->second);
            }
        }

        /*!
         * \brief Get the communicated sparsity pattern
         * \return the vector with the sparsity pattern
         */
        std::vector<std::set<int> >& sparsity ()
        {
            return sparsity_;
        }

    private:
        const DofMapper& mapper_;
        const std::map<IdType,int>& gid2Index_;
        const std::map<int,IdType>& index2GID_;
        std::vector<std::set<int> > sparsity_;
        Matrix& A_;
        const ParallelISTLHelper<GridView, LinearSolverTraits>& helper_;

    }; // class MatPatternExchange

    //! Local matrix blocks associated with the global id set
    struct MatEntry
    {
        IdType first;
        BlockType second;
        MatEntry (const IdType& f, const BlockType& s) : first(f), second(s) {}
        MatEntry () {}
    };

    //! A DataHandle class to exchange matrix entries
    class MatEntryExchange
        : public Dune::CommDataHandleIF<MatEntryExchange,MatEntry>
    {
        using RowIterator = typename Matrix::RowIterator;
        using ColIterator = typename Matrix::ColIterator;
    public:
        //! Export type of data for message buffer
        using DataType = MatEntry;

        /*!
         * \brief Constructor
         * \param[in] mapper The local dof mapper.
         * \param[in] g2i Global to local index map.
         * \param[in] i2g Local to global index map.
         * \param[in] A Matrix to operate on.
         */
        MatEntryExchange (const DofMapper& mapper, const std::map<IdType,int>& g2i,
                          const std::map<int,IdType>& i2g,
                          Matrix& A)
            : mapper_(mapper), gid2Index_(g2i), index2GID_(i2g), A_(A)
        {}

        /*!
         * \brief Returns true if data for given valid codim should be communicated
         */
        bool contains (int dim, int codim) const
        {
            return (codim==dofCodim);
        }

        /*!
         * \brief Returns true if size of data per entity of given dim and codim is a constant
         */
        bool fixedsize (int dim, int codim) const
        {
            return false;
        }

        /*!
         * \brief How many objects of type DataType have to be sent for a given entity
         */
        template<class EntityType>
        size_t size (EntityType& e) const
        {
            int i = mapper_.index(e);
            int n = 0;
            for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
                typename std::map<int,IdType>::const_iterator it = index2GID_.find(j.index());
                if (it != index2GID_.end())
                    n++;
            }

            return n;
        }

        /*!
         * \brief Pack data from user to message buffer
         */
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
            int i = mapper_.index(e);
            for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
                typename std::map<int,IdType>::const_iterator it=index2GID_.find(j.index());
                if (it != index2GID_.end())
                    buff.write(MatEntry(it->second,*j));
            }

        }

        /*!
         * \brief Unpack data from message buffer to user
         */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
            int i = mapper_.index(e);
            for (size_t k = 0; k < n; k++)
            {
                MatEntry m;
                buff.read(m);
                // only add entries corresponding to border entities
                typename std::map<IdType,int>::const_iterator it = gid2Index_.find(m.first);
                if (it != gid2Index_.end())
                    if (A_[i].find(it->second) != A_[i].end())
                        A_[i][it->second] += m.second;
            }
        }

    private:
        const DofMapper& mapper_;
        const std::map<IdType,int>& gid2Index_;
        const std::map<int,IdType>& index2GID_;
        Matrix& A_;

    }; // class MatEntryExchange

    [[deprecated("Use extendMatrix instead. Will be removed after 3.1!")]]
    void getExtendedMatrix (Matrix& A, const ParallelISTLHelper<GridView, LinearSolverTraits>& helper)
    { extendMatrix(A, helper); }

    /*!
     * \brief communicates values for the sparsity pattern of the new matrix.
     * \param A Matrix to operate on.
     * \param helper ParallelelISTLHelper.
     */
    void extendMatrix(Matrix& A, const ParallelISTLHelper<GridView, LinearSolverTraits>& helper)
    {
        if (gridView_.comm().size() > 1)
        {
            Matrix tmp(A);
            std::size_t nnz = 0;
            // get entries from other processes
            MatPatternExchange datahandle(mapper_, gid2Index_, index2GID_, A, helper);
            gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface,
                                              Dune::ForwardCommunication);
            std::vector<std::set<int> >& sparsity = datahandle.sparsity();
            // add own entries, count number of nonzeros
            for (RowIterator i = A.begin(); i != A.end(); ++i)
            {
                for (ColIterator j = A[i.index()].begin(); j != A[i.index()].end(); ++j)
                {
                    if (sparsity[i.index()].find(j.index()) == sparsity[i.index()].end())
                        sparsity[i.index()].insert(j.index());
                }
                nnz += sparsity[i.index()].size();
            }

            A.setSize(tmp.N(), tmp.N(), nnz);
            A.setBuildMode(Matrix::row_wise);
            typename Matrix::CreateIterator citer = A.createbegin();
            using Iter = typename std::vector<std::set<int> >::const_iterator;
            for (Iter i = sparsity.begin(), end = sparsity.end(); i!=end; ++i, ++citer)
            {
                using SIter = typename std::set<int>::const_iterator;
                for (SIter si = i->begin(), send = i->end(); si!=send; ++si)
                    citer.insert(*si);
            }

            // set matrix old values
            A = 0;
            for (RowIterator i = tmp.begin(); i != tmp.end(); ++i)
                for (ColIterator j = tmp[i.index()].begin(); j != tmp[i.index()].end(); ++j)
                    A[i.index()][j.index()] = tmp[i.index()][j.index()];
        }
    }

    /*!
     * \brief Sums up the entries corresponding to border vertices.
     * \param A Matrix to operate on.
     */
    void sumEntries (Matrix& A)
    {
        if (gridView_.comm().size() > 1)
        {
            MatEntryExchange datahandle(mapper_, gid2Index_, index2GID_, A);
            gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface,
                                              Dune::ForwardCommunication);
        }
    }

#if HAVE_MPI
    /*!
     * \brief Extends the sparsity pattern of the discretization matrix for AMG.
     * \param A A reference to the matrix to change.
     */
    [[deprecated("Use extendMatrix instead. Will be removed after 3.1!")]]
    void getExtendedMatrix (Matrix& A)
    { extendMatrix(A); }
    /*!
     * \brief Extends the sparsity pattern of the discretization matrix for AMG.
     * \param A A reference to the matrix to change.
     */
    void extendMatrix(Matrix& A)
    {
        if (gridView_.comm().size() > 1)
        {
            Matrix tmp(A);
            std::size_t nnz = 0;
            // get entries from other processes
            MatPatternExchange datahandle(mapper_, gid2Index_, index2GID_, A, *this);
            gridView_.communicate(datahandle, Dune::InteriorBorder_InteriorBorder_Interface,
                                  Dune::ForwardCommunication);
            std::vector<std::set<int> >& sparsity = datahandle.sparsity();

            // add own entries, count number of nonzeros
            for (RowIterator i = A.begin(); i != A.end(); ++i)
            {
                for (ColIterator j = A[i.index()].begin(); j != A[i.index()].end(); ++j)
                {
                    if (sparsity[i.index()].find(j.index()) == sparsity[i.index()].end())
                        sparsity[i.index()].insert(j.index());
                }
                nnz += sparsity[i.index()].size();
            }

            A.setSize(tmp.N(), tmp.N(), nnz);
            A.setBuildMode(Matrix::row_wise);
            typename Matrix::CreateIterator citer = A.createbegin();
            using Iter = typename std::vector<std::set<int> >::const_iterator;
            for (Iter i = sparsity.begin(), end = sparsity.end(); i!=end; ++i, ++citer){
                using SIter = typename std::set<int>::const_iterator;
                for (SIter si = i->begin(), send = i->end(); si!=send; ++si)
                    citer.insert(*si);
            }

            // set matrix old values
            A = 0;
            for (RowIterator i = tmp.begin(); i != tmp.end(); ++i)
                for (ColIterator j = tmp[i.index()].begin(); j != tmp[i.index()].end(); ++j)
                    A[i.index()][j.index()] = tmp[i.index()][j.index()];

            sumEntries(A);
        }
    }
#endif

private:
    const GridView gridView_;
    const DofMapper& mapper_;
    std::map<IdType, int> gid2Index_;
    std::map<int, IdType> index2GID_;

}; // class EntityExchanger

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
template<class GridView, class LinearSolverTraits, bool isParallel>
struct LinearAlgebraPreparator
{
    using DofMapper = typename LinearSolverTraits::DofMapper;
    using ParallelHelper = ParallelISTLHelper<GridView, LinearSolverTraits>;
    using Comm = typename LinearSolverTraits::Comm;
    using LinearOperator = typename LinearSolverTraits::LinearOperator;
    using ScalarProduct = typename LinearSolverTraits::ScalarProduct;

    template<class Matrix, class Vector>
    static void prepareLinearAlgebra(Matrix& A, Vector& b,
                                     int& rank,
                                     std::shared_ptr<Comm>& comm,
                                     std::shared_ptr<LinearOperator>& fop,
                                     std::shared_ptr<ScalarProduct>& sp,
                                     ParallelHelper& pHelper,
                                     const bool firstCall)
    {
        comm = std::make_shared<Comm>();
        fop = std::make_shared<LinearOperator>(A);
        sp = std::make_shared<ScalarProduct>();
    }
};

#if HAVE_MPI
/*!
 * \brief Specialization for the parallel case.
 */
template<class GridView, class LinearSolverTraits>
struct LinearAlgebraPreparator<GridView, LinearSolverTraits, true>
{
    using DofMapper = typename LinearSolverTraits::DofMapper;
    using ParallelHelper = ParallelISTLHelper<GridView, LinearSolverTraits>;
    using Comm = typename LinearSolverTraits::Comm;
    using LinearOperator = typename LinearSolverTraits::LinearOperator;
    using ScalarProduct = typename LinearSolverTraits::ScalarProduct;

    template<class Matrix, class Vector>
    static void prepareLinearAlgebra(Matrix& A, Vector& b,
                                     int& rank,
                                     std::shared_ptr<Comm>& comm,
                                     std::shared_ptr<LinearOperator>& fop,
                                     std::shared_ptr<ScalarProduct>& sp,
                                     ParallelHelper& pHelper,
                                     const bool firstCall)
    {
        const auto category = LinearSolverTraits::isNonOverlapping ?
                              Dune::SolverCategory::nonoverlapping : Dune::SolverCategory::overlapping;

        if (LinearSolverTraits::isNonOverlapping && firstCall)
            pHelper.initGhostsAndOwners();

        comm = std::make_shared<Comm>(pHelper.gridView().comm(), category);

        if (LinearSolverTraits::isNonOverlapping)
        {
            // extend the matrix pattern such that it is usable for the parallel solver
            EntityExchanger<GridView, LinearSolverTraits> exchanger(pHelper.gridView(), pHelper.dofMapper());
            exchanger.extendMatrix(A, pHelper);
            exchanger.sumEntries(A);
        }
        pHelper.createParallelIndexSet(*comm);

        fop = std::make_shared<LinearOperator>(A, *comm);
        sp = std::make_shared<ScalarProduct>(*comm);
        rank = comm->communicator().rank();

        // Make rhs consistent
        if (LinearSolverTraits::isNonOverlapping)
            pHelper.makeNonOverlappingConsistent(b);
    }
};

#endif // HAVE_MPI

} // end namespace Dumux

#endif // DUMUX_PARALLELHELPERS_HH
