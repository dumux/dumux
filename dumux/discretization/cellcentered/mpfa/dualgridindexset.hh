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

#include <dune/common/reservedvector.hh>

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
 * \tparam maxE The maximum admissible number of elements around vertices.
 * \tparam maxB The maximum admissible number of branches on intersections.
 *              This is only to be specified for network grids and defaults to 1
 *              for normal grids.
 */
template< class GI, class LI, int dim, int maxE, int maxB = 2 >
class CCMpfaDualGridNodalIndexSet
{
    using DimIndexVector = Dune::ReservedVector<GI, dim>;

public:
    //! Export the used index types
    using GridIndexType = GI;
    using LocalIndexType = LI;

    //! Export the specified maximum admissible sizes
    static constexpr int dimension = dim;
    static constexpr int maxBranches = maxB;
    static constexpr int maxNumElementsAtNode = maxE*(maxBranches-1);
    static constexpr int maxNumScvfsAtNode = maxNumElementsAtNode*dim;

    //! Export the stencil types used
    using GridStencilType = Dune::ReservedVector<GridIndexType, maxNumElementsAtNode>;
    using LocalStencilType = Dune::ReservedVector<LocalIndexType, maxNumElementsAtNode>;

    //! Export the type used for storing the global scvf indices at this node
    using GridScvfStencilType = Dune::ReservedVector<GridIndexType, maxNumScvfsAtNode>;

    //! Data structure to store the neighboring scv indices of an scvf (grid/local indices)
    using ScvfNeighborIndexSet = Dune::ReservedVector<GridIndexType, maxBranches>;
    using ScvfNeighborLocalIndexSet = Dune::ReservedVector<LocalIndexType, maxBranches>;

    //! Constructor
    CCMpfaDualGridNodalIndexSet() : numBoundaryScvfs_(0) {}

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
        ScvfNeighborIndexSet scvIndices;
        scvIndices.push_back(insideScvIdx);
        for (const auto& outsideIdx : outsideScvIndices)
                scvIndices.push_back(outsideIdx);

        // if scvf is on boundary, increase counter
        if (boundary)
            numBoundaryScvfs_++;

        // insert data on the new scv
        scvfIndices_.push_back(scvfIdx);
        scvfIsOnBoundary_.push_back(boundary);
        scvfNeighborScvIndices_.push_back(scvIndices);

        // if entry for the inside scv exists, append scvf local index, create otherwise
        auto it = std::find( scvIndices_.begin(), scvIndices_.end(), insideScvIdx );
        if (it != scvIndices_.end())
            localScvfIndicesInScv_[ std::distance(scvIndices_.begin(), it) ].push_back(curScvfLocalIdx);
        else
        {
            localScvfIndicesInScv_.push_back({});
            localScvfIndicesInScv_.back().push_back(curScvfLocalIdx);
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
    const GridStencilType& globalScvIndices() const { return scvIndices_; }

    //! returns the global scvf indices connected to this dual grid node
    const GridScvfStencilType& globalScvfIndices() const { return scvfIndices_; }

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
    const ScvfNeighborIndexSet& neighboringScvIndices(unsigned int i) const
    { return scvfNeighborScvIndices_[i]; }

private:
    GridStencilType scvIndices_;                                                       //!< The indices of the scvs around a dual grid node
    Dune::ReservedVector<DimIndexVector, maxNumElementsAtNode> localScvfIndicesInScv_; //!< Maps to each scv a list of scvf indices embedded in it

    GridScvfStencilType scvfIndices_;                             //!< the indices of the scvfs around a dual grid node
    std::size_t numBoundaryScvfs_;                                //!< stores how many boundary scvfs are embedded in this dual grid node
    Dune::ReservedVector<bool, maxNumScvfsAtNode> scvfIsOnBoundary_; //!< Maps to each scvf a boolean to indicate if it is on the boundary

    //! The neighboring scvs for the scvfs (order: 0 - inside, 1..n - outside)
    Dune::ReservedVector<ScvfNeighborIndexSet, maxNumScvfsAtNode> scvfNeighborScvIndices_;
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
    CCMpfaDualGridIndexSet(const GridView& gridView) : nodalIndexSets_(gridView.size(NodalIndexSet::dimension)) {}

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
