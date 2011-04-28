// $Id$
/*****************************************************************************
 *   Copyright (C) 2011 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A BCRS matrix which creates an algebraic overlap of
 *        arbitrary size.
 */
#ifndef DUMUX_OVERLAP_MATRIX_HH
#define DUMUX_OVERLAP_MATRIX_HH

#include <dune/grid/common/datahandleif.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <list>
#include <set>


namespace Dumux {

/*!
 * \brief Uses communication on the grid to find the initial seed list
 *        of indices.
 *
 * \todo implement this class generically. For this, it must be
 *       possible to query the mapper whether it contains entities of
 *       a given codimension without the need to hand it an actual
 *       entity.
 */
template <class GridView, class VertexMapper>
class VertexSeedListFromGrid : public Dune::CommDataHandleIF< VertexSeedListFromGrid<GridView, VertexMapper>,
                                                              int >
{
    typedef int ProcessRank;
    typedef int RowIndex;
    typedef std::pair<RowIndex, ProcessRank> RowProcessPair;
    typedef std::list<RowProcessPair> SeedList;

public:
    VertexSeedListFromGrid(const GridView &gv, 
                           const VertexMapper &map)
        : gv_(gv), map_(map)
    {
        gv.communicate(*this, 
                       Dune::InteriorBorder_InteriorBorder_Interface, 
                       Dune::ForwardCommunication);
    };

    // data handle methods
    bool contains (int dim, int codim) const
    { return dim == codim; }
    
    bool fixedsize(int dim, int codim) const
    { return true; }
    
    template<class EntityType> 
    size_t size(const EntityType &e) const
    { return 1; };

    template<class MessageBufferImp, class EntityType> 
    void gather(MessageBufferImp &buff, const EntityType &e) const 
    { buff.write(gv_.comm().rank()); };

    template<class MessageBufferImp, class EntityType> 
    void scatter(MessageBufferImp &buff, const EntityType &e, size_t n)
    {
        int remoteRank;
        buff.read(remoteRank);
        int vertIdx = map_.map(e);
        seedList_.push_back(RowProcessPair(vertIdx, remoteRank));
    };
    
    const SeedList &seedList() const
    { return seedList_; }

private:
    const GridView gv_;
    const VertexMapper &map_;
    SeedList seedList_;
};


/*!
 * \brief A BCRS matrix which creates an algebraic overlap of
 *        arbitrary size.
 *
 * This class can be used for the parallel linear solver.
 */
template<class BlockType>
class OverlapMatrix : public Dune::BCRSMatrix<BlockType>
{
    typedef Dune::BCRSMatrix<BlockType> ParentType;
    typedef typename BlockType::field_type Scalar;

    OverlapMatrix(const OverlapMatrix &A)
    {}

public:
    typedef int ProcessRank;
    typedef int RowIndex;
    typedef std::pair<RowIndex, ProcessRank> RowProcessPair;
    typedef std::list<RowProcessPair> SeedList;

    OverlapMatrix(const ParentType &A,
                  const SeedList &seedList,
                  int overlapSize = 1)
    {
        overlapSet_.resize(A.N());
        
        appendOverlapIndices_(A, seedList, overlapSize);
        createMatrix_(A);
    };

    void printOverlap()
    {
        int n = overlapSet_.size();
        for (int i = 0; i < n; ++i) {
            const std::set<int> &rowSet = overlapSet_[i];
            if (rowSet.size() == 0)
                continue;
            
            std::cout << "Overlap processes for row " << i << ": ";
            
            std::set<int>::const_iterator it = rowSet.begin();
            std::set<int>::const_iterator endIt = rowSet.end();
            for (; it != endIt; ++it) {
                std::cout << *it << " ";
            };
            std::cout << "\n";
        };
    };
                  
protected:
    void appendOverlapIndices_(const ParentType &A,
                               const SeedList &seedList,
                               int overlapSize)
    {
        // add all processes in the seed rows of the current overlap
        // level
        SeedList::const_iterator it = seedList.begin();
        SeedList::const_iterator endIt = seedList.end();
        for (; it != endIt; ++it)
            overlapSet_[it->first].insert(it->second);

        // if we have an overlap of 0, we're finished and break the
        // recursion
        if (overlapSize == 0)
            return;

        // find the seed list for the next overlap level using the
        // seed set for the current level
        SeedList nextSeedList;
        it = seedList.begin();
        for (; it != endIt; ++it) {
            int rowIdx = it->first;
            int processIdx = it->second;

            // find all column indices in the row. The indices of the
            // columns are the additional indices of the overlap which
            // we would like to add
            typedef typename ParentType::ConstColIterator ColIterator;
            ColIterator colIt = A[rowIdx].begin();
            ColIterator colEndIt = A[rowIdx].end();
            for (; colIt != colEndIt; ++colIt) {
                int newIdx = colIt.index();

                // if the process is already is in the overlap of the
                // column index, ignore this column index!
                if (overlapSet_[newIdx].count(processIdx) > 0)
                    continue;
                
                // add the current processes to the seed list for the
                // next overlap level
                nextSeedList.push_back(RowProcessPair(newIdx, processIdx));
            }
        }                          
        
        // Perform the same excercise for the next overlap level
        appendOverlapIndices_(A, nextSeedList, overlapSize - 1);
    };
    
    void createMatrix_(const ParentType &A)
    {
    };

    // stores the set of process ranks which are in the overlap for a
    // given row index.
    std::vector<std::set<int> > overlapSet_;
};

} // namespace Dumux

#endif
