// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 *
 * \brief Provides a helper class for nonoverlapping
 *        decomposition using the PDELab AMG.
 */
#ifndef DUMUX_AMGPARALLELHELPERS_HH
#define DUMUX_AMGPARALLELHELPERS_HH

#include <dune/common/version.hh>

#include <dumux/implicit/box/boxproperties.hh>
#include <dumux/implicit/cellcentered/ccproperties.hh>
#include <dumux/decoupled/common/pressureproperties.hh>
#include <dumux/linear/amgproperties.hh>

namespace Dumux
{

/*!
 * \brief A parallel helper class providing a nonoverlapping
 *        decomposition of all degrees of freedom
 */
// operator that resets result to zero at constrained DOFS
template<class TypeTag>
class ParallelISTLHelper
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP(TypeTag, AMGLocalFemMap) LocalFemMap;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum { dofCodim = LocalFemMap::dofCodim };

    class BaseGatherScatter
    {
    public:
        BaseGatherScatter(const Problem& problem)
        : problem_(problem)
        {}

        template<class EntityType>
        int map(const EntityType& e) const
        {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            return problem_.model().dofMapper().index(e);
#else
            return problem_.model().dofMapper().map(e);
#endif
        }

    private:
        const Problem& problem_;
    };

    /**
     * @brief GatherScatter implementation that makes a right hand side in the box model consistent.
     */
    template<class V>
    class ConsistencyBoxGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<ConsistencyBoxGatherScatter<V>,typename V::block_type>
    {
    public:
        typedef typename V::block_type DataType;
        
        ConsistencyBoxGatherScatter(V& container, const Problem& problem)
            : BaseGatherScatter(problem), container_(container)
        {}
        
        bool contains(int dim, int codim) const
        {
            return dofCodim==codim;
        }

        bool fixedsize(int dim, int codim) const
        {
            return true;
        }

        template<class EntityType>
        size_t size (EntityType& e) const
        {
            return 1;
        }

        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
            buff.write(container_[this->map(e)]);
        }

        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
            typename V::block_type block;
            buff.read(block);
            container_[this->map(e)]+=block;
        }
    private:
        V& container_;
    };
    

    /**
     * @brief Writes 1<<24 to each data item (of the container) that is gathered or scattered
     * and is neither interior nor border.
     *
     * Can be used to mark ghost cells.
     */
    class GhostGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<GhostGatherScatter,std::size_t>
    {
    public:
        typedef std::size_t DataType;

        GhostGatherScatter(std::vector<std::size_t>& ranks,
                           const Problem& problem)
            : BaseGatherScatter(problem), ranks_(ranks)
        {}


        bool contains(int dim, int codim) const
        {
            return dofCodim==codim;
        }

        bool fixedsize(int dim, int codim) const
        {
            return true;
        }

        template<class EntityType>
        size_t size (EntityType& e) const
        {
            return 1;
        }

        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
            std::size_t& data= ranks_[this->map(e)];
            if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
                data = (1<<24);
            buff.write(data);
        }

        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
            std::size_t x;
            std::size_t& data = ranks_[this->map(e)];
            buff.read(x);
            if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
                data= (1<<24);
        }
    private:
        std::vector<std::size_t>& ranks_;
    };

    /**
     * @brief GatherScatter handle that sets 1<<24 for data items neither associated to
     * the interior or border and take the minimum when scattering.
     *
     * Used to compute an owner rank for each unknown.
     */
    class InteriorBorderGatherScatter
        : public BaseGatherScatter, 
          public Dune::CommDataHandleIF<InteriorBorderGatherScatter,std::size_t>
    {
    public:
        typedef std::size_t DataType;

        InteriorBorderGatherScatter(std::vector<std::size_t>& ranks,
                                    const Problem& problem)
            : BaseGatherScatter(problem), ranks_(ranks)
        {}


        bool contains(int dim, int codim) const
        {
            return dofCodim==codim;
        }

        bool fixedsize(int dim, int codim) const
        {
            return true;

        }

        template<class EntityType>
        size_t size (EntityType& e) const
        {
            return 1;
        }

        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {

            std::size_t& data = ranks_[this->map(e)];
            if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
                data = (1<<24);
            buff.write(data);
        }

        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
            std::size_t x;
            std::size_t& data = ranks_[this->map(e)];
            buff.read(x);
            if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
                data = x;
            else
                data = std::min(data,x);
        }
    private:
        std::vector<std::size_t>& ranks_;
    };

    /**
     * @brief GatherScatter handle for finding out about neighbouring processor ranks.
     *
     */
    struct NeighbourGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<NeighbourGatherScatter,int>
    {
        typedef int DataType;

        NeighbourGatherScatter(const Problem& problem, int rank_, std::set<int>& neighbours_)
            : BaseGatherScatter(problem), myrank(rank_), neighbours(neighbours_)
        {}


        bool contains(int dim, int codim) const
        {
            return dofCodim==codim;
        }

        bool fixedsize(int dim, int codim) const
        {
            return true;
        }

        template<class EntityType>
        size_t size (EntityType& e) const
        {
            return 1;
        }

        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType &e) const
        {
            buff.write(myrank);
        }

        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType &e, size_t n)
        {
            int x;
            buff.read(x);
            neighbours.insert(x);
        }
        int myrank;
        std::set<int>& neighbours;
    };


    /**
     * @brief GatherScatter handle for finding out about neighbouring processor ranks.
     *
     */
    struct SharedGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<SharedGatherScatter,int>
    {
        typedef int DataType;

        SharedGatherScatter(std::vector<int>& shared,
                            const Problem& problem)
            : BaseGatherScatter(problem), shared_(shared)
        {}

        bool contains(int dim, int codim) const
        {
            return dofCodim==codim;
        }

        bool fixedsize(int dim, int codim) const
        {
            return true;

        }

        template<class EntityType>
        size_t size (EntityType& e) const
        {
            return 1;
        }

        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, EntityType& e) const
        {
            int data=true;
            buff.write(data);
        }

        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType &e, size_t n)
        {
            int x;
            buff.read(x);
            int& data= shared_[this->map(e)];
            data = data || x;
        }
    private:
        std::vector<int>& shared_;

    };

    /**
     * @brief GatherScatter handle for finding out about neighbouring processor ranks.
     *
     */
    template<typename GI>
    struct GlobalIndexGatherScatter
        : public BaseGatherScatter,
          public Dune::CommDataHandleIF<GlobalIndexGatherScatter<GI>, GI>
    {
        typedef GI DataType;
        GlobalIndexGatherScatter(std::vector<GI>& gindices,
                                 const Problem& problem)
            : BaseGatherScatter(problem), gindices_(gindices)
        {}

        bool contains(int dim, int codim) const
        {
            return dofCodim==codim;
        }

        bool fixedsize(int dim, int codim) const
        {
            return true;
        }

        template<class EntityType>
        size_t size (EntityType& e) const
        {
            return 1;
        }

        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
            buff.write(gindices_[this->map(e)]);
        }

        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
            DataType x;
            buff.read(x);
            gindices_[this->map(e)] = std::min(gindices_[this->map(e)], x);
        }
    private:
        std::vector<GI>& gindices_;
    };

public:

    ParallelISTLHelper (const Problem& problem, int verbose=1)
        : problem_(problem), verbose_(verbose), initialized_(false)
    {}
    
    // \brief Initializes the markers for ghosts and owners with the correct size and values.
    //
    void initGhostsAndOwners(){
        owner_.resize(problem_.model().dofMapper().size(),
                      problem_.gridView().comm().rank());
        isGhost_.resize(problem_.model().dofMapper().size(),0.0);
        // find out about ghosts
        GhostGatherScatter ggs(owner_,problem_);

        if (problem_.gridView().comm().size()>1)
            problem_.gridView().communicate(ggs,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);

        isGhost_ = owner_;

        // partition interior/border
        InteriorBorderGatherScatter dh(owner_, problem_);

        if (problem_.gridView().comm().size()>1)
            problem_.gridView().communicate(dh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);

        // convert vector into mask vector
        for(auto v=owner_.begin(), vend=owner_.end(); v!=vend;++v)
            if(*v==problem_.gridView().comm().rank())
                *v=1.0;
            else
                *v=0.0;

        initialized_=true;
    }

    // keep only DOFs assigned to this processor
    template<typename W>
    void mask (W& w) const
    {
        auto v1=w.begin();

        for(auto v2=owner_.begin(), vend=owner_.end(); v2!=vend;++v1,++v2)
            v1*=v2;
    }

    // access to mask vector
    double mask (std::size_t i) const
    {
        return owner_[i];
    }

    // access to ghost vector
    double ghost (std::size_t i) const
    {
        return isGhost_[i];
    }
    
    // \brief Make a vector of the box model consistent.
    template<typename B, typename A>
    void makeNonOverlappingConsistent(Dune::BlockVector<B,A>& v)
    {
#if HAVE_MPI
        const GridView& gridview = problem_.gridView();
        ConsistencyBoxGatherScatter<Dune::BlockVector<B,A> > gs(v, problem_);
        if (gridview.comm().size()>1)
            gridview.communicate(gs,Dune::InteriorBorder_InteriorBorder_Interface,
                                 Dune::ForwardCommunication);
#endif
    }
    

#if HAVE_MPI

    /**
     * @brief Creates a matrix suitable for parallel AMG and the parallel information
     *
     *
     * @tparam MatrixType The type of the ISTL matrix used.
     * @tparam Comm The type of the OwnerOverlapCopyCommunication
     * @param m The local matrix.
     * @param c The parallel information object providing index set, interfaces and
     * communicators.
     */
    template<typename MatrixType, typename Comm>
    void createIndexSetAndProjectForAMG(MatrixType& m, Comm& c);
#endif
private:
    const Problem& problem_;
    std::vector<std::size_t> owner_; // vector to identify unique decomposition
    std::vector<std::size_t> isGhost_; //vector to identify ghost dofs
    int verbose_; //verbosity
    // \brief whether isGhost and owner arrays are initialized
    bool initialized_;
};

/**
 * @brief Helper class for adding up matrix entries on border.
 * @tparam GridOperator The grid operator to work on.
 * @tparam MatrixType The MatrixType.
 */
template<class TypeTag>
class EntityExchanger
{
    typedef typename GET_PROP_TYPE(TypeTag, JacobianMatrix) JacobianMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP(TypeTag, AMGLocalFemMap) LocalFemMap;
    enum { numEq = JacobianMatrix::block_type::rows};
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar,numEq,numEq> > Matrix;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum {dim = GridView::dimension};
    typedef typename GridView::Traits::Grid Grid;
    typedef typename Matrix::block_type BlockType;
    typedef typename GridView::template Codim<LocalFemMap::dofCodim>::Iterator EntityIterator;
    typedef typename Grid::Traits::GlobalIdSet IDS;
    typedef typename IDS::IdType IdType;
    typedef typename Matrix::RowIterator RowIterator;
    typedef typename Matrix::ColIterator ColIterator;

    enum { dofCodim = LocalFemMap::dofCodim };

public:
    /*! \brief Constructor. Sets up the local to global relations.
      \param[in] problem The problem to be solved.
    */
    EntityExchanger(const Problem& problem)
    : problem_(problem)
    {
        gid2Index_.clear();
        index2GID_.clear();

        const GridView& gridView = problem.gridView();

        EntityIterator entityEndIt = gridView.template end<LocalFemMap::dofCodim>();
        for (EntityIterator entityIt = gridView.template begin<LocalFemMap::dofCodim>();
             entityIt != entityEndIt; ++entityIt)
        {
            if (entityIt->partitionType() == Dune::BorderEntity)
            {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
                int localIdx = problem_.model().dofMapper().index(*entityIt);
#else
                int localIdx = problem_.model().dofMapper().map(*entityIt);
#endif
                IdType dofIdxGlobal = gridView.grid().globalIdSet().id(*entityIt);

                std::pair<IdType,int> g2iPair(dofIdxGlobal, localIdx);
                gid2Index_.insert(g2iPair);

                std::pair<int,IdType> i2gPair(localIdx, dofIdxGlobal);
                index2GID_.insert(i2gPair);

            }
        }
    }

    /**
     * @brief A DataHandle class to exchange matrix sparsity patterns.
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
        : public Dune::CommDataHandleIF<MatPatternExchange,IdType> {
        typedef typename Matrix::RowIterator RowIterator;
        typedef typename Matrix::ColIterator ColIterator;
    public:
        //! Export type of data for message buffer
        typedef IdType DataType;

        /** @brief Returns true if data for given valid codim should be communicated
         */
        bool contains (int dim, int codim) const
        {
            return (codim==dofCodim);
        }

        /** @brief Returns true if size of data per entity of given dim and codim is a constant
         */
        bool fixedsize (int dim, int codim) const
        {
            return false;
        }

        /** @brief How many objects of type DataType have to be sent for a given entity
         */
        template<class EntityType>
        size_t size (EntityType& e) const
        {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            int i = problem_.model().dofMapper().index(e);
#else
            int i = problem_.model().dofMapper().map(e);
#endif
            int n = 0;
            for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
                typename std::map<int,IdType>::const_iterator it = index2GID_.find(j.index());
                if (it != index2GID_.end())
                    n++;
            }

            return n;
        }

        /** @brief Pack data from user to message buffer
         */
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            int i = problem_.model().dofMapper().index(e);
#else
            int i = problem_.model().dofMapper().map(e);
#endif
            for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
                typename std::map<int,IdType>::const_iterator it=index2GID_.find(j.index());
                if (it != index2GID_.end())
                    buff.write(it->second);
            }

        }

        /** @brief Unpack data from message buffer to user
         */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            int i = problem_.model().dofMapper().index(e);
#else
            int i = problem_.model().dofMapper().map(e);
#endif
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

        /**
         * @brief Get the communicated sparsity pattern
         * @return the vector with the sparsity pattern
         */
        std::vector<std::set<int> >& sparsity ()
        {
            return sparsity_;
        }

        /** @brief Constructor
            @param[in] problem The problem to be solved.
            @param[in] g2i Global to local index map.
            @param[in] i2g Local to global index map.
            @param[in] A Matrix to operate on.
            @param[in] helper parallel istl helper.
        */
        MatPatternExchange (const Problem& problem,
                            const std::map<IdType,int>& g2i,
                            const std::map<int,IdType>& i2g, Matrix& A,
                            const ParallelISTLHelper<TypeTag>& helper)
            : problem_(problem), gid2Index_(g2i), index2GID_(i2g),
              sparsity_(A.N()), A_(A), helper_(helper)
        {}

    private:
        const Problem& problem_;
        const std::map<IdType,int>& gid2Index_;
        const std::map<int,IdType>& index2GID_;
        std::vector<std::set<int> > sparsity_;
        Matrix& A_;
        const ParallelISTLHelper<TypeTag>& helper_;
    };

    //! Local matrix blocks associated with the global id set
    struct MatEntry
    {
        IdType first;
        BlockType second;
        MatEntry (const IdType& f, const BlockType& s) : first(f),second(s) {}
        MatEntry () {}
    };

    //! A DataHandle class to exchange matrix entries
    class MatEntryExchange
        : public Dune::CommDataHandleIF<MatEntryExchange,MatEntry> {
        typedef typename Matrix::RowIterator RowIterator;
        typedef typename Matrix::ColIterator ColIterator;
    public:
        //! Export type of data for message buffer
        typedef MatEntry DataType;

        /** @brief Returns true if data for given valid codim should be communicated
         */
        bool contains (int dim, int codim) const
        {
            return (codim==dofCodim);
        }

        /** @brief Returns true if size of data per entity of given dim and codim is a constant
         */
        bool fixedsize (int dim, int codim) const
        {
            return false;
        }

        /** @brief How many objects of type DataType have to be sent for a given entity
         */
        template<class EntityType>
        size_t size (EntityType& e) const
        {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            int i = problem_.model().dofMapper().index(e);
#else
            int i = problem_.model().dofMapper().map(e);
#endif
            int n = 0;
            for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
                typename std::map<int,IdType>::const_iterator it = index2GID_.find(j.index());
                if (it != index2GID_.end())
                    n++;
            }

            return n;
        }

        /** @brief Pack data from user to message buffer
         */
        template<class MessageBuffer, class EntityType>
        void gather (MessageBuffer& buff, const EntityType& e) const
        {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            int i = problem_.model().dofMapper().index(e);
#else
            int i = problem_.model().dofMapper().map(e);
#endif
            for (ColIterator j = A_[i].begin(); j != A_[i].end(); ++j)
            {
                typename std::map<int,IdType>::const_iterator it=index2GID_.find(j.index());
                if (it != index2GID_.end())
                    buff.write(MatEntry(it->second,*j));
            }

        }

        /** @brief Unpack data from message buffer to user
         */
        template<class MessageBuffer, class EntityType>
        void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
        {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
            int i = problem_.model().dofMapper().index(e);
#else
            int i = problem_.model().dofMapper().map(e);
#endif
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

        /** @brief Constructor
            @param[in] problem The problem to be solved.
            @param[in] g2i Global to local index map.
            @param[in] i2g Local to global index map.
            @param[in] A Matrix to operate on.
        */
        MatEntryExchange (const Problem& problem, const std::map<IdType,int>& g2i,
                          const std::map<int,IdType>& i2g,
                          Matrix& A)
            : problem_(problem), gid2Index_(g2i), index2GID_(i2g), A_(A)
        {}

    private:
        const Problem& problem_;
        const std::map<IdType,int>& gid2Index_;
        const std::map<int,IdType>& index2GID_;
        Matrix& A_;
    };

    /** @brief communicates values for the sparsity pattern of the new matrix.
        @param A Matrix to operate on.
        @param helper ParallelelISTLHelper.
    */
    void getExtendedMatrix (Matrix& A,const ParallelISTLHelper<TypeTag>& helper)
    {
        if (problem_.gridView().comm().size() > 1) {
            Matrix tmp(A);
            std::size_t nnz=0;
            // get entries from other processes
            MatPatternExchange datahandle(problem_, gid2Index_, index2GID_, A, helper);
            problem_.gridView().communicate(datahandle,
                                                    Dune::InteriorBorder_InteriorBorder_Interface,
                                                    Dune::ForwardCommunication);
            std::vector<std::set<int> >& sparsity = datahandle.sparsity();
            // add own entries, count number of nonzeros
            for (RowIterator i = A.begin(); i != A.end(); ++i){
                for (ColIterator j = A[i.index()].begin(); j != A[i.index()].end(); ++j){
                    if (sparsity[i.index()].find(j.index()) == sparsity[i.index()].end())
                        sparsity[i.index()].insert(j.index());
                }
                nnz += sparsity[i.index()].size();
            }
            A.setSize(tmp.N(), tmp.N(), nnz);
            A.setBuildMode(Matrix::row_wise);
            typename Matrix::CreateIterator citer = A.createbegin();
            typedef typename std::vector<std::set<int> >::const_iterator Iter;
            for (Iter i = sparsity.begin(), end = sparsity.end(); i!=end; ++i, ++citer){
                typedef typename std::set<int>::const_iterator SIter;
                for (SIter si = i->begin(), send = i->end(); si!=send; ++si)
                    citer.insert(*si);
            }
            // set matrix old values
            A = 0;
            for (RowIterator i = tmp.begin(); i != tmp.end(); ++i)
                for (ColIterator j = tmp[i.index()].begin(); j != tmp[i.index()].end(); ++j){
                    A[i.index()][j.index()] = tmp[i.index()][j.index()];
                }
        }
    }

    /** @brief Sums up the entries corresponding to border vertices.
        @param A Matrix to operate on.
    */
    void sumEntries (Matrix& A)
    {
        if (problem_.gridView().comm().size() > 1)
        {
            MatEntryExchange datahandle(problem_, gid2Index_, index2GID_, A);
            problem_.gridView().communicate(datahandle,
                                                    Dune::InteriorBorder_InteriorBorder_Interface,
                                                    Dune::ForwardCommunication);
        }
    }

#if HAVE_MPI
    /**
     * @brief Extends the sparsity pattern of the discretization matrix for AMG.
     * @param A A reference to the matrix to change.
     */
    void getExtendedMatrix (Matrix& A) const;
#endif

private:
    const Problem& problem_;
    std::map<IdType,int> gid2Index_;
    std::map<int,IdType> index2GID_;
};

#if HAVE_MPI
template<class TypeTag>
void EntityExchanger<TypeTag>::getExtendedMatrix (Matrix& A) const
{
    const GridView& gridView = problem_.gridView();
    if (gridView.comm().size() > 1) {
        Matrix tmp(A);
        std::size_t nnz=0;
        // get entries from other processes
        MatPatternExchange datahandle(problem_, gid2Index_, index2GID_, A, *this);
        gridView.communicate(datahandle,
                             Dune::InteriorBorder_InteriorBorder_Interface,
                             Dune::ForwardCommunication);
        std::vector<std::set<int> >& sparsity = datahandle.sparsity();
        // add own entries, count number of nonzeros
        for (RowIterator i = A.begin(); i != A.end(); ++i){
            for (ColIterator j = A[i.index()].begin(); j != A[i.index()].end(); ++j){
                if (sparsity[i.index()].find(j.index()) == sparsity[i.index()].end())
                    sparsity[i.index()].insert(j.index());
            }
            nnz += sparsity[i.index()].size();
        }
        A.setSize(tmp.N(), tmp.N(), nnz);
        A.setBuildMode(Matrix::row_wise);
        typename Matrix::CreateIterator citer = A.createbegin();
        typedef typename std::vector<std::set<int> >::const_iterator Iter;
        for (Iter i = sparsity.begin(), end = sparsity.end(); i!=end; ++i, ++citer){
            typedef typename std::set<int>::const_iterator SIter;
            for (SIter si = i->begin(), send = i->end(); si!=send; ++si)
                citer.insert(*si);
        }
        // set matrix old values
        A = 0;
        for (RowIterator i = tmp.begin(); i != tmp.end(); ++i)
            for (ColIterator j = tmp[i.index()].begin(); j != tmp[i.index()].end(); ++j){
                A[i.index()][j.index()] = tmp[i.index()][j.index()];
            }
        sumEntries(A);
    }
}

template<class TypeTag>
template<typename M, typename C>
void ParallelISTLHelper<TypeTag>::createIndexSetAndProjectForAMG(M& m, C& c)
{
    if(!initialized_){
        // This is the first time this function is called.
        // Therefore we need to initialize the marker vectors for ghosts and
        // owned dofs
        initGhostsAndOwners();
    }
    const GridView& gridview = problem_.gridView();

    // First find out which dofs we share with other processors
    std::vector<int> sharedDofs(problem_.model().dofMapper().size(), false);

    SharedGatherScatter sgs(sharedDofs, problem_);

    if (gridview.comm().size()>1)
        gridview.communicate(sgs,Dune::All_All_Interface,
                             Dune::ForwardCommunication);

    // Count shared dofs that we own
    typedef typename C::ParallelIndexSet::GlobalIndex GlobalIndex;
    GlobalIndex count=0;
    auto owned=owner_.begin();
    
    for(auto v=sharedDofs.begin(), vend=sharedDofs.end(); v != vend; ++v, ++owned)
        if(*v && *owned==1.0)
            ++count;

    Dune::dverb<<gridview.comm().rank()<<": shared count is "<< count.touint()
         <<std::endl;

    std::vector<GlobalIndex> counts(gridview.comm().size());
    gridview.comm().allgather(&count, 1, &(counts[0]));

    // Compute start index start_p = \sum_{i=0}^{i<p} counts_i
    GlobalIndex start=0;
    for(int i=0; i<gridview.comm().rank(); ++i)
        start=start+counts[i];
    //std::cout<<gv.comm().rank()<<": start index = "<<start.touint()<<std::endl;


    typedef std::vector<GlobalIndex> GIVector;
    GIVector scalarIndices(problem_.model().dofMapper().size(),
                           std::numeric_limits<GlobalIndex>::max());

    auto shared=sharedDofs.begin();
    auto index=scalarIndices.begin();

    for(auto i=owner_.begin(), iend=owner_.end(); i!=iend; ++i, ++shared, ++index)
        if(*i==1.0 && *shared){
            *index=start;
            ++start;
        }

    // publish global indices for the shared DOFS to other processors.
    typedef GlobalIndexGatherScatter<GlobalIndex> GIGS;
    GIGS gigs(scalarIndices, problem_);
    if (gridview.comm().size()>1)
        gridview.communicate(gigs,Dune::All_All_Interface,
                             Dune::ForwardCommunication);


    // Setup the index set
    c.indexSet().beginResize();
    index=scalarIndices.begin();
    auto ghost=isGhost_.begin();

    for(auto i=owner_.begin(), iend=owner_.end(); i!=iend; ++i, ++ghost, ++index)
    {
        Dune::OwnerOverlapCopyAttributeSet::AttributeSet attr;
        if(*index!=std::numeric_limits<GlobalIndex>::max()){
            // global index exist in index set
            if(*i>0){
                // This dof is managed by us.
                attr = Dune::OwnerOverlapCopyAttributeSet::owner;
            }
            else if ( *ghost==(1<<24) && ( c.getSolverCategory() ==
                                           static_cast<int>(Dune::SolverCategory::nonoverlapping)) ){
                //use attribute overlap for ghosts in novlp grids
                attr = Dune::OwnerOverlapCopyAttributeSet::overlap;
            }
            else {
                attr = Dune::OwnerOverlapCopyAttributeSet::copy;
            }
            c.indexSet().add(*index, typename C::ParallelIndexSet::LocalIndex(i-owner_.begin(), attr));
        }
    }
    c.indexSet().endResize();
    //std::cout<<gv.comm().rank()<<": index set size = "<<c.indexSet().size()<<std::endl;
    //std::cout<<gv.comm().rank()<<": "<<c.indexSet()<<std::endl;

    // Compute neighbours using communication
    typedef NeighbourGatherScatter NeighbourGS;
    std::set<int> neighbours;
    NeighbourGatherScatter ngs(problem_, gridview.comm().rank(),
                               neighbours);

    if (gridview.comm().size()>1)
        gridview.communicate(ngs,Dune::All_All_Interface,
                             Dune::ForwardCommunication);
    c.remoteIndices().setNeighbours(neighbours);
    //std::cout<<gv.comm().rank()<<": no neighbours="<<neighbours.size()<<std::endl;

    c.remoteIndices().template rebuild<false>();
    //std::cout<<c.remoteIndices()<<std::endl;

}
#endif
} // end namespace Dumux

#endif
