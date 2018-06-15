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
 * \brief copydoc Dumux::FacetGridVertexAdapter
 */
#ifndef DUMUX_FACETCOUPLING_FACET_GRID_VERTEX_ADAPTER_HH
#define DUMUX_FACETCOUPLING_FACET_GRID_VERTEX_ADAPTER_HH

#include <dune/grid/common/mcmgmapper.hh>

namespace Dumux {

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Adapter to be plugged on top of a grid creator to allow obtaining
 *        vertex index within a d-dimensional grid (corresponding to the given
 *        bulk grid id) for a given (d-1)-dimensional grid vertex.
 *
 * \tparam GridCreator A grid creator containing a hierarchy of face-conforming grids.
 * \tparam bulkId The grid id of the d-dimensional grid within the hierarchy
 * \tparam facetId The grid if of the (d-1)-dimensional grid within the hierarchy
 */
template<class GridCreator, int bulkId = 0, int facetId = 1>
class FacetGridVertexAdapter
{
    using FacetGrid = typename GridCreator::template Grid<facetId>;
    using FacetGridView = typename FacetGrid::LeafGridView;
    using FacetGridVertex = typename FacetGridView::template Codim<FacetGrid::dimension>::Entity;
    using FacetGridElement = typename FacetGridView::template Codim<0>::Entity;
    using FacetGridIndexType = typename FacetGridView::IndexSet::IndexType;

    using BulkGrid = typename GridCreator::template Grid<bulkId>;
    using BulkGridView = typename BulkGrid::LeafGridView;
    using BulkGridVertex = typename BulkGridView::template Codim<BulkGrid::dimension>::Entity;
    using BulkIndexType = typename BulkGridView::IndexSet::IndexType;

    // Find out if the given bulk grid is the one with highest dimensionality among the created grids
    static constexpr bool bulkHasHighestDimension = (int(BulkGrid::dimension) == GridCreator::bulkDim);

    // check if provided id combination makes sense
    static_assert( int(FacetGrid::dimension) == int(BulkGrid::dimension) - 1,
                   "Grid dimension mismatch! Please check the provided domain ids!" );
    static_assert( int(FacetGrid::dimensionworld) == int(BulkGrid::dimensionworld),
                   "Grid world dimension mismatch! All grids must have the same world dimension" );

public:

    //! The constructor
    FacetGridVertexAdapter(const GridCreator& gridCreator)
    : gridCreatorPtr_(&gridCreator)
    {
        using BulkMapper = Dune::MultipleCodimMultipleGeomTypeMapper<BulkGridView>;
        const auto& bulkGridView = gridCreator.template grid<bulkId>().leafGridView();
        auto bulkVertexMapper = BulkMapper(bulkGridView, Dune::mcmgVertexLayout());

        // insertion to grid index map
        bulkInsertionToGridVIdx_.resize(bulkGridView.size(BulkGrid::dimension));
        const auto& bulkGridFactory = gridCreatorPtr_->template gridFactory<bulkId>();
        for (const auto& v : vertices(bulkGridView))
            bulkInsertionToGridVIdx_[ bulkGridFactory.insertionIndex(v) ] = bulkVertexMapper.index(v);

        // set up index map for the bulk vertex indices
        makeBulkIndexMap_(gridCreator);

        // determine which bulk vertices lie on facet elements
        bulkVertexIsOnFacetGrid_.resize(bulkGridView.size(BulkGrid::dimension), false);
        const auto& facetGridView = gridCreator.template grid<facetId>().leafGridView();
        const auto& facetGridFactory = gridCreatorPtr_->template gridFactory<facetId>();
        for (const auto& v : vertices(facetGridView))
        {
            const auto insIdx = facetGridFactory.insertionIndex(v);
            const auto highestLevelInsIdx = gridCreatorPtr_->lowDimVertexIndices(facetId)[insIdx];
            const auto bulkIdx = bulkHasHighestDimension
                                 ? bulkInsertionToGridVIdx_[ highestLevelInsIdx ]
                                 : bulkInsertionToGridVIdx_[ highestLevelInsertionToBulkInsertionIdx_[ highestLevelInsIdx] ];
            bulkVertexIsOnFacetGrid_[ bulkIdx ] = true;
        }
    }

    /*!
     * \brief Returns the index within the d-dimensional grid of a vertex
     *        of the (d-1)-dimensional grid.
     * \note  Leads to undefined behaviour if called for a vertex which doens't
     *        exist on the d-dimensional grid
     */
    BulkIndexType bulkGridVertexIndex(const FacetGridVertex& v) const
    {
        const auto insIdx = gridCreatorPtr_->template gridFactory<facetId>().insertionIndex(v);
        const auto highestLevelInsIdx = gridCreatorPtr_->lowDimVertexIndices(facetId)[insIdx];

        if (bulkHasHighestDimension)
            return bulkInsertionToGridVIdx_[ highestLevelInsIdx ];
        else
            return bulkInsertionToGridVIdx_[ highestLevelInsertionToBulkInsertionIdx_[ highestLevelInsIdx] ];
    }

    /*!
     * \brief Returns true if the vertex of the d-dimensional grid with
     *        the given vertex index also exists on the (d-1)-dimensional grid
     */
    bool isOnFacetGrid(const BulkGridVertex& v) const
    {
        const auto bulkInsIdx = gridCreatorPtr_->template gridFactory<bulkId>().insertionIndex(v);
        const auto bulkVIdx = bulkHasHighestDimension
                              ? bulkInsertionToGridVIdx_[ bulkInsIdx ]
                              : bulkInsertionToGridVIdx_[ highestLevelInsertionToBulkInsertionIdx_[ bulkInsIdx] ];
        return bulkVertexIsOnFacetGrid_[ bulkVIdx ];
    }

    /*!
     * \brief Returns true if the vertex of the d-dimensional grid with
     *        the given vertex index also exists on the (d-1)-dimensional grid
     */
    bool isOnFacetGrid(BulkIndexType bulkVIdx) const
    { return bulkVertexIsOnFacetGrid_[ bulkVIdx ]; }

    /*!
     * \brief Returns true if a (d-1)-dimensional element is embedded in
     *        the d-dimensional domain
     */
    bool isEmbedded(const FacetGridElement& e) const
    { return numEmbedments(e) > 0; }

    /*!
     * \brief Returns the number of d-dimensional elements in which the
     *        given (d-1)-dimensional element is embedded in
     */
    std::size_t numEmbedments(const FacetGridElement& e) const
    { return gridCreatorPtr_->template embedmentEntityIndices<facetId>(e).size(); }

private:
    //! Determine the map from the insertion idx of the highest-dimensional grid to bulk insertion index
    template< bool isHighest = bulkHasHighestDimension, std::enable_if_t<!isHighest, int> = 0 >
    void makeBulkIndexMap_(const GridCreator& gridCreator)
    {
        const auto& highLevelGridView = gridCreator.template grid<0>().leafGridView();
        highestLevelInsertionToBulkInsertionIdx_.resize(highLevelGridView.size(GridCreator::bulkDim));

        const auto& bulkVertexIndices = gridCreatorPtr_->lowDimVertexIndices(bulkId);
        for (std::size_t bulkInsIdx = 0; bulkInsIdx < bulkVertexIndices.size(); ++bulkInsIdx)
            highestLevelInsertionToBulkInsertionIdx_[ bulkVertexIndices[bulkInsIdx] ] = bulkInsIdx;
    }

    //! If the given bulk grid is on the highest level of grid creation, we do not need the map
    template< bool isHighest = bulkHasHighestDimension, std::enable_if_t<isHighest, int> = 0 >
    void makeBulkIndexMap_(const GridCreator& gridCreator) {}

    // data members
    const GridCreator* gridCreatorPtr_;
    std::vector<bool> bulkVertexIsOnFacetGrid_;
    std::vector<BulkIndexType> bulkInsertionToGridVIdx_;
    std::vector<BulkIndexType> highestLevelInsertionToBulkInsertionIdx_;
};

} // end namespace Dumux

#endif
