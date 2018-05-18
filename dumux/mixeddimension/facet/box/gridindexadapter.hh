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
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief copydoc Dumux::FacetCouplingVertexMapper
 */
#ifndef DUMUX_FACETCOUPLING_BOX_FACET_GRIDINDEX_ADAPTER_HH
#define DUMUX_FACETCOUPLING_BOX_FACET_GRIDINDEX_ADAPTER_HH

#include <dune/grid/common/mcmgmapper.hh>

namespace Dumux
{

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Adapter to be plugged on top of a grid creator to allow obtaining
 *        the index of a vertex within  a d-dimensional grid for a given
 *        (d-1)-dimensional grid vertex.
 *
 * \tparam GridCreator A grid creator containing a hierarchy of face-conforming grids.
 * \tparam bulkId The grid id of the d-dimensional grid within the hierarchy
 * \tparam lowDimId The grid if of the (d-1)-dimensional grid within the hierarchy
 */
template<class GridCreator, int bulkId, int lowDimId>
class FacetGridIndexAdapter
{
    using LowDimGrid = typename GridCreator::template Grid<lowDimId>;
    using LowDimGridView = typename LowDimGrid::LeafGridView;
    using LowDimGridVertex = typename LowDimGridView::template Codim<LowDimGrid::dimension>::Entity;
    using LowDimGridElement = typename LowDimGridView::template Codim<0>::Entity;
    using LowDimIndexType = typename LowDimGridView::IndexSet::IndexType;

    using BulkGrid = typename GridCreator::template Grid<bulkId>;
    using BulkGridView = typename BulkGrid::LeafGridView;
    using BulkIndexType = typename BulkGridView::IndexSet::IndexType;

    // Find out if the given bulk grid is on the highest level among the created grids
    static constexpr bool isHighestLevel = (int(BulkGrid::dimension) == GridCreator::bulkDim);

    // check if provided id combination makes sense
    static_assert( int(LowDimGrid::dimension) == int(BulkGrid::dimension) - 1,
                   "Grid dimension mismatch! Please check the provided domain ids!" );
    static_assert( int(LowDimGrid::dimensionworld) == int(BulkGrid::dimensionworld),
                   "Grid world dimension mismatch! All grids must have the same world dimension" );

public:

    //! the constructor
    FacetGridIndexAdapter(const GridCreator& gridCreator)
    : gridCreatorPtr_(&gridCreator)
    {
        using BulkMapper = Dune::MultipleCodimMultipleGeomTypeMapper<BulkGridView>;
        const auto& bulkGridView = gridCreator.template grid<bulkId>().leafGridView();
        const auto& bulkGridFactory = gridCreatorPtr_->template gridFactory<bulkId>();
        auto bulkVertexMapper = BulkMapper(bulkGridView, Dune::mcmgVertexLayout());

        // insertion to grid index map
        bulkInsertionToGridVIdx_.resize(bulkGridView.size(BulkGrid::dimension));
        for (const auto& v : vertices(bulkGridView))
            bulkInsertionToGridVIdx_[ bulkGridFactory.insertionIndex(v) ] = bulkVertexMapper.index(v);

        // set up index map for the bulk vertex indices
        makeBulkIndexMap_(gridCreator);

        // determine which bulk vertices lie on facet elements
        bulkVertexIsOnLowDimGrid_.resize(bulkGridView.size(BulkGrid::dimension), false);
        const auto& lowDimGridView = gridCreator.template grid<lowDimId>().leafGridView();
        const auto& lowDimGridFactory = gridCreatorPtr_->template gridFactory<lowDimId>();
        for (const auto& v : vertices(lowDimGridView))
        {
            const auto insIdx = lowDimGridFactory.insertionIndex(v);
            const auto highestLevelInsIdx = gridCreatorPtr_->lowDimVertexIndices(lowDimId)[insIdx];
            const auto bulkInsIdx = isHighestLevel
                                    ? bulkInsertionToGridVIdx_[ highestLevelInsIdx ]
                                    : bulkInsertionToGridVIdx_[ highestLevelInsertionToBulkInsertionIdx_[ highestLevelInsIdx] ];
            bulkVertexIsOnLowDimGrid_[ bulkInsIdx ] = true;
        }
    }

    //! returns the index within the d-dimensional grid of a vertex of the (d-1)-dimensional grid
    //! leads to undefined behaviour if this called for a vertex which doesn't exist on the d-dimensional grid
    BulkIndexType index(const LowDimGridVertex& v) const
    {
        const auto insIdx = gridCreatorPtr_->template gridFactory<lowDimId>().insertionIndex(v);
        const auto highestLevelInsIdx = gridCreatorPtr_->lowDimVertexIndices(lowDimId)[insIdx];

        if (isHighestLevel)
            return bulkInsertionToGridVIdx_[ highestLevelInsIdx ];
        else
            return bulkInsertionToGridVIdx_[ highestLevelInsertionToBulkInsertionIdx_[ highestLevelInsIdx] ];
    }

    //! returns true if a (d-1)-dimensional element is embedded in the d-dimensional domain
    bool isEmbedded(const LowDimGridElement& e) const
    { return gridCreatorPtr_->template embedmentEntityIndices<lowDimId>(e).size() > 0; }

    //! returns true if a bulk grid vertex lies on a lowDim grid element
    bool vertexIsOnLowDimGrid(BulkIndexType bulkVIdx) const
    { return bulkVertexIsOnLowDimGrid_[ bulkVIdx ]; }

private:
    //! Determine the map from the insertion idx of the highest-dimensional grid to bulk insertion index
    template< bool isHighest = isHighestLevel, std::enable_if_t<!isHighest, int> = 0 >
    void makeBulkIndexMap_(const GridCreator& gridCreator)
    {
        const auto& highLevelGridView = gridCreator.template grid<0>().leafGridView();
        highestLevelInsertionToBulkInsertionIdx_.resize(highLevelGridView.size(GridCreator::bulkDim));

        const auto& bulkVertexIndices = gridCreatorPtr_->lowDimVertexIndices(bulkId);
        for (std::size_t bulkInsIdx = 0; bulkInsIdx < bulkVertexIndices.size(); ++bulkInsIdx)
            highestLevelInsertionToBulkInsertionIdx_[ bulkVertexIndices[bulkInsIdx] ] = bulkInsIdx;
    }

    //! If the given bulk grid is on the highest level of grid creation, we do not need the map
    template< bool isHighest = isHighestLevel, std::enable_if_t<isHighest, int> = 0 >
    void makeBulkIndexMap_(const GridCreator& gridCreator) {}

    // data members
    const GridCreator* gridCreatorPtr_;
    std::vector<bool> bulkVertexIsOnLowDimGrid_;
    std::vector<BulkIndexType> bulkInsertionToGridVIdx_;
    std::vector<BulkIndexType> highestLevelInsertionToBulkInsertionIdx_;
};
} // end namespace Dumux

#endif
