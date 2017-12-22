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
 * \ingroup CCMpfaDiscretization
 * \brief Class for the index sets of the dual grid in mpfa schemes.
 */
#ifndef DUMUX_DISCRETIZATION_MPFA_DUALGRID_INDEX_SET_HH
#define DUMUX_DISCRETIZATION_MPFA_DUALGRID_INDEX_SET_HH

#include <cassert>
#include <vector>
#include <algorithm>

namespace Dumux
{

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Nodal index set for mpfa schemes, constructed
 *        around grid vertices.
 *
 * \tparam GI The type used for indices on the grid
 * \tparam LI The type used for indexing in interaction volumes
 * \tparam dim The dimension of the grid
 */
template< class GI, class LI, int dim >
class DualGridNodalIndexSet
{
public:
    using GridIndexType = GI;
    using LocalIndexType = LI;

    // we use dynamic containers to store indices here
    using GridIndexContainer = std::vector<GridIndexType>;
    using LocalIndexContainer = std::vector<LocalIndexType>;

    //! Constructor
    DualGridNodalIndexSet() : numBoundaryScvfs_(0) {}

    //! Inserts data for a given scvf
    template<typename SubControlVolumeFace>
    void insert(const SubControlVolumeFace& scvf)
    {
        insert(scvf.boundary(),
               scvf.index(),
               scvf.insideScvIdx(),
               scvf.outsideScvIndices());
    }

    //! Inserts scvf data
    void insert(const bool boundary,
                const GridIndexType scvfIdx,
                const GridIndexType insideScvIdx,
                const std::vector<GridIndexType>& outsideScvIndices)
    {
        // this should always be called only once per scvf
        assert(std::find(scvfIndices_.begin(), scvfIndices_.end(), scvfIdx ) == scvfIndices_.end()
               && "scvf has already been inserted!");

        // the local index of the scvf data about to be inserted
        const LocalIndexType curScvfLocalIdx = scvfIndices_.size();

        // add grid scvf data
        GridIndexContainer scvIndices;
        scvIndices.reserve(outsideScvIndices.size()+1);
        scvIndices.push_back(insideScvIdx);
        scvIndices.insert(scvIndices.end(), outsideScvIndices.begin(), outsideScvIndices.end());

        // if scvf is on boundary, increase counter
        if (boundary)
            numBoundaryScvfs_++;

        // insert data on the new scv
        scvfIndices_.push_back(scvfIdx);
        scvfIsOnBoundary_.push_back(boundary);
        scvfNeighborScvIndices_.emplace_back( std::move(scvIndices) );

        // if entry for the inside scv exists, append scvf local index, create otherwise
        auto it = std::find( scvIndices_.begin(), scvIndices_.end(), insideScvIdx );
        if (it != scvIndices_.end())
            localScvfIndicesInScv_[ std::distance(scvIndices_.begin(), it) ].push_back(curScvfLocalIdx);
        else
        {
            LocalIndexContainer localScvfIndices;
            localScvfIndices.reserve(dim);
            localScvfIndices.push_back(curScvfLocalIdx);
            localScvfIndicesInScv_.emplace_back( std::move(localScvfIndices) );
            scvIndices_.push_back(insideScvIdx);
        }
    }

    //! returns the number of scvs around the node
    std::size_t numScvs() const { return scvIndices_.size(); }

    //! returns the number of scvfs around the node
    std::size_t numScvfs() const { return scvfIndices_.size(); }

    //! returns the number of boundary scvfs around the node
    std::size_t numBoundaryScvfs() const { return numBoundaryScvfs_; }

    //! returns the global scv indices connected to this dual grid node
    const GridIndexContainer& globalScvIndices() const { return scvIndices_; }

    //! returns the global scvf indices connected to this dual grid node
    const GridIndexContainer& globalScvfIndices() const { return scvfIndices_; }

    //! returns whether or not the i-th scvf is on a domain boundary
    bool scvfIsOnBoundary(unsigned int i) const { return scvfIsOnBoundary_[i]; }

    //! returns the global scv idx of the i-th scv
    GridIndexType scvIdxGlobal(unsigned int i) const { return scvIndices_[i]; }

    //! returns the index of the i-th scvf
    GridIndexType scvfIdxGlobal(unsigned int i) const { return scvfIndices_[i]; }

    //! returns the global index of the j-th scvf embedded in the i-th scv
    GridIndexType scvfIdxGlobal(unsigned int i, unsigned int j) const
    {
        assert(j < dim); // only dim faces can be embedded in an scv
        return scvfIndices_[ localScvfIndicesInScv_[i][j] ];
    }

    //! returns the node-local index of the j-th scvf embedded in the i-th scv
    LocalIndexType scvfIdxLocal(unsigned int i, unsigned int j) const
    {
        assert(j < dim); // only dim faces can be embedded in an scv
        return localScvfIndicesInScv_[i][j];
    }

    //! returns the indices of the neighboring scvs of the i-th scvf
    const GridIndexContainer& neighboringScvIndices(unsigned int i) const
    { return scvfNeighborScvIndices_[i]; }

private:
    GridIndexContainer scvIndices_;                          //!< The indices of the scvs around a dual grid node
    std::vector<LocalIndexContainer> localScvfIndicesInScv_; //!< Maps to each scv a list of scvf indices embedded in it

    GridIndexContainer scvfIndices_;                         //!< the indices of the scvfs around a dual grid node
    std::vector<bool> scvfIsOnBoundary_;                     //!< Maps to each scvf a boolean to indicate if it is on the boundary
    std::vector<GridIndexContainer> scvfNeighborScvIndices_; //!< The neighboring scvs for the scvfs (order: 0 - inside, 1..n - outside)
    std::size_t numBoundaryScvfs_;                           //!< stores how many boundary scvfs are embedded in this dual grid node
};

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Class for the index sets of the dual grid in mpfa schemes.
 *
 * \tparam GridIndexType The type used for indices on the grid
 * \tparam LocalIndexType The type used for indexing in interaction volumes
 * \tparam dim The dimension of the grid
 */
template< class GridIndexType, class LocalIndexType, int dim >
class CCMpfaDualGridIndexSet
{
public:
    using NodalIndexSet = DualGridNodalIndexSet< GridIndexType, LocalIndexType, dim >;

    //! Default constructor should not be used
    CCMpfaDualGridIndexSet() = delete;

    //! Constructor taking a grid view
    template< class GridView >
    CCMpfaDualGridIndexSet(const GridView& gridView) : nodalIndexSets_(gridView.size(dim)) {}

    //! Access with an scvf
    template< class SubControlVolumeFace >
    const NodalIndexSet& operator[] (const SubControlVolumeFace& scvf) const
    { return nodalIndexSets_[scvf.vertexIndex()]; }

    template< class SubControlVolumeFace >
    NodalIndexSet& operator[] (const SubControlVolumeFace& scvf)
    { return nodalIndexSets_[scvf.vertexIndex()]; }

    //! Access with an index
    const NodalIndexSet& operator[] (GridIndexType i) const { return nodalIndexSets_[i]; }
    NodalIndexSet& operator[] (GridIndexType i) { return nodalIndexSets_[i]; }

private:
    std::vector<NodalIndexSet> nodalIndexSets_;
};

} // end namespace Dumux

#endif
