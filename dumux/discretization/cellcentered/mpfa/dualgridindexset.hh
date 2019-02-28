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
 * \ingroup CCMpfaDiscretization
 * \brief Class for the index sets of the dual grid in mpfa schemes.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_DUALGRID_INDEX_SET_HH
#define DUMUX_DISCRETIZATION_MPFA_DUALGRID_INDEX_SET_HH

#include <cassert>
#include <vector>
#include <algorithm>

#include <dune/common/reservedvector.hh>
#include <dumux/common/indextraits.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Default traits to be used in conjuntion
 *        with the dual grid nodal index set.
 *
 * \tparam GV The grid view type
 */
template<class GV>
struct NodalIndexSetDefaultTraits
{
    using GridView = GV;
    using GridIndexType = typename IndexTraits<GV>::GridIndex;
    using LocalIndexType = typename IndexTraits<GV>::LocalIndex;

    //! per default, we use dynamic data containers (iv size unknown)
    template< class T > using NodalScvDataStorage = std::vector< T >;
    template< class T > using NodalScvfDataStorage = std::vector< T >;

    //! store data on neighbors of scvfs in static containers if possible
    template< class T >
    using ScvfNeighborDataStorage = typename std::conditional_t< (int(GV::dimension)<int(GV::dimensionworld)),
                                                                 std::vector< T >,
                                                                 Dune::ReservedVector< T, 2 > >;
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Nodal index set for mpfa schemes, constructed
 *        around grid vertices.
 *
 * \tparam T The traits class to be used
 */
template< class T >
class CCMpfaDualGridNodalIndexSet
{
    using LI = typename T::LocalIndexType;
    using GI = typename T::GridIndexType;

    using DimIndexVector = Dune::ReservedVector<LI, T::GridView::dimension>;
    using ScvfIndicesInScvStorage = typename T::template NodalScvDataStorage< DimIndexVector >;

public:
    //! Export the traits type
    using Traits = T;

    //! Export the index types used
    using LocalIndexType = LI;
    using GridIndexType = GI;

    //! Export the stencil types used
    using NodalGridStencilType = typename T::template NodalScvDataStorage< GI >;
    using NodalLocalStencilType = typename T::template NodalScvDataStorage< LI >;
    using NodalGridScvfStencilType = typename T::template NodalScvfDataStorage< GI >;

    //! Data structure to store the neighboring scv indices of an scvf (grid/local indices)
    using ScvfNeighborLocalIndexSet = typename T::template ScvfNeighborDataStorage< LI >;

    //! Constructor
    CCMpfaDualGridNodalIndexSet() : numBoundaryScvfs_(0) {}

    //! Inserts data for a given scvf
    template<typename SubControlVolumeFace>
    void insert(const SubControlVolumeFace& scvf)
    {
        insert(scvf.index(), scvf.insideScvIdx(), scvf.boundary());
    }

    //! Inserts scvf data
    void insert(const GridIndexType scvfIdx,
                const GridIndexType insideScvIdx,
                const bool boundary)
    {
        // this should always be called only once per scvf
        assert(std::count(scvfIndices_.begin(), scvfIndices_.end(), scvfIdx ) == 0 && "scvf already inserted!");

        // the local index of the scvf data about to be inserted
        const LocalIndexType curScvfLocalIdx = scvfIndices_.size();

        // if scvf is on boundary, increase counter
        if (boundary) numBoundaryScvfs_++;

        // insert data on the new scv
        scvfIndices_.push_back(scvfIdx);
        scvfIsOnBoundary_.push_back(boundary);

        // if entry for the inside scv exists append data, create otherwise
        auto it = std::find( scvIndices_.begin(), scvIndices_.end(), insideScvIdx );
        if (it != scvIndices_.end())
        {
            const auto scvIdxLocal = std::distance(scvIndices_.begin(), it);
            scvfInsideScvIndices_.push_back(scvIdxLocal);
            localScvfIndicesInScv_[scvIdxLocal].push_back(curScvfLocalIdx);
        }
        else
        {
            scvfInsideScvIndices_.push_back(scvIndices_.size());
            localScvfIndicesInScv_.push_back({curScvfLocalIdx});
            scvIndices_.push_back(insideScvIdx);
        }
    }

    //! returns the number of scvs around the node
    std::size_t numScvs() const
    { return scvIndices_.size(); }

    //! returns the number of scvfs around the node
    std::size_t numScvfs() const
    { return scvfIndices_.size(); }

    //! returns the number of boundary scvfs around the node
    std::size_t numBoundaryScvfs() const
    { return numBoundaryScvfs_; }

    //! returns the grid scv indices connected to this dual grid node
    const NodalGridStencilType& gridScvIndices() const
    { return scvIndices_; }

    //! returns the grid scvf indices connected to this dual grid node
    const NodalGridScvfStencilType& gridScvfIndices() const
    { return scvfIndices_; }

    //! returns whether or not the i-th scvf is on a domain boundary
    bool scvfIsOnBoundary(unsigned int i) const
    {
        assert(i < numScvfs());
        return scvfIsOnBoundary_[i];
    }

    //! returns the grid scv idx of the i-th scv
    GridIndexType gridScvIndex(unsigned int i) const
    {
        assert(i < numScvs());
        return scvIndices_[i];
    }

    //! returns the index of the i-th scvf
    GridIndexType gridScvfIndex(unsigned int i) const
    {
        assert(i < numScvfs());
        return scvfIndices_[i];
    }

    //! returns the grid index of the j-th scvf embedded in the i-th scv
    GridIndexType gridScvfIndex(unsigned int i, unsigned int j) const
    {
        assert(i < numScvs());
        assert(j < localScvfIndicesInScv_[i].size());
        return scvfIndices_[ localScvfIndicesInScv_[i][j] ];
    }

    //! returns the node-local index of the j-th scvf embedded in the i-th scv
    LocalIndexType localScvfIndex(unsigned int i, unsigned int j) const
    {
        assert(i < numScvs());
        assert(j < localScvfIndicesInScv_[i].size());
        return localScvfIndicesInScv_[i][j];
    }

    //! returns the node-local index of the inside scv of the i-th scvf
    LocalIndexType insideScvLocalIndex(unsigned int i) const
    {
        assert(i < numScvfs());
        return scvfInsideScvIndices_[i];
    }

private:
    NodalGridStencilType scvIndices_;               //!< The indices of the scvs around a dual grid node
    ScvfIndicesInScvStorage localScvfIndicesInScv_; //!< Maps to each scv a list of scvf indices embedded in it

    std::size_t numBoundaryScvfs_;                                         //!< stores how many boundary scvfs are embedded in this dual grid node
    NodalGridScvfStencilType scvfIndices_;                                 //!< the indices of the scvfs around a dual grid node
    typename T::template NodalScvfDataStorage< bool > scvfIsOnBoundary_;   //!< Maps to each scvf a boolean to indicate if it is on the boundary
    typename T::template NodalScvfDataStorage< LI > scvfInsideScvIndices_; //!< The inside local scv index for each scvf
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class for the index sets of the dual grid in mpfa schemes.
 *
 * \tparam NI The type used for the nodal index sets.
 */
template< class NI >
class CCMpfaDualGridIndexSet
{
public:
    using NodalIndexSet = NI;
    using GridIndexType = typename NodalIndexSet::GridIndexType;

    //! Default constructor should not be used
    CCMpfaDualGridIndexSet() = delete;

    //! Constructor taking a grid view
    template< class GridView >
    CCMpfaDualGridIndexSet(const GridView& gridView)
    : nodalIndexSets_(gridView.size(GridView::dimension))
    {}

    //! Access with an scvf
    template< class SubControlVolumeFace >
    const NodalIndexSet& operator[] (const SubControlVolumeFace& scvf) const
    { return nodalIndexSets_[scvf.vertexIndex()]; }

    template< class SubControlVolumeFace >
    NodalIndexSet& operator[] (const SubControlVolumeFace& scvf)
    { return nodalIndexSets_[scvf.vertexIndex()]; }

    //! Access with an index
    const NodalIndexSet& operator[] (GridIndexType i) const
    { return nodalIndexSets_[i]; }

    NodalIndexSet& operator[] (GridIndexType i)
    { return nodalIndexSets_[i]; }

private:
    std::vector<NodalIndexSet> nodalIndexSets_;
};

} // end namespace Dumux

#endif
