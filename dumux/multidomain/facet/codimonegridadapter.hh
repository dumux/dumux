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
 * \ingroup FacetCoupling
 * \brief copydoc Dumux::CodimOneGridAdapter
 */
#ifndef DUMUX_FACETCOUPLING_CODIM_ONE_GRID_ADAPTER_HH
#define DUMUX_FACETCOUPLING_CODIM_ONE_GRID_ADAPTER_HH

#include <memory>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/indextraits.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Adapter that  allows retrieving information on a d-dimensional grid for
 *        entities of a (d-1)-dimensional grid. This lower-dimensional grid is
 *        assumed to be facet-conforming to the d-dimensional grid. This class
 *        can be used in the context of models where a sub-domain lives on the
 *        facets of a bulk grid.
 *
 * \tparam Embeddings Class containing the embedments of entites between grids of codim 1
 * \tparam bulkGridId The grid id of the d-dimensional grid within the hierarchy
 * \tparam facetGridId The grid id of the (d-1)-dimensional grid within the hierarchy
 */
template<class Embeddings, int bulkGridId = 0, int facetGridId = 1>
class CodimOneGridAdapter
{
    // Extract some types of the facet-conforming grid of codimension one
    using FacetGridView = typename Embeddings::template GridView<facetGridId>;
    using FacetGridVertex = typename FacetGridView::template Codim<FacetGridView::dimension>::Entity;
    using FacetGridElement = typename FacetGridView::template Codim<0>::Entity;
    using FacetGridIndexType = typename IndexTraits<FacetGridView>::GridIndex;

    // Extract some types of the bulk grid
    using BulkGridView = typename Embeddings::template GridView<bulkGridId>;
    using BulkMapper = Dune::MultipleCodimMultipleGeomTypeMapper<BulkGridView>;
    using BulkReferenceElements = typename Dune::ReferenceElements<typename BulkGridView::ctype, BulkGridView::dimension>;
    using BulkGridElement = typename BulkGridView::template Codim<0>::Entity;
    using BulkGridIntersection = typename BulkGridView::Intersection;
    using BulkGridVertex = typename BulkGridView::template Codim<BulkGridView::dimension>::Entity;
    using BulkIndexType = typename IndexTraits<BulkGridView>::GridIndex;

    // Find out if the given bulk grid is the one with highest dimensionality among the created grids
    static constexpr bool bulkHasHighestDimension = (int(BulkGridView::dimension) == Embeddings::bulkDim);

    // check if provided id combination makes sense
    static_assert( int(FacetGridView::dimension) == int(BulkGridView::dimension) - 1,
                   "Grid dimension mismatch! Please check the provided domain ids!" );
    static_assert( int(FacetGridView::dimensionworld) == int(BulkGridView::dimensionworld),
                   "Grid world dimension mismatch! All grids must have the same world dimension" );

public:

    //! The constructor
    CodimOneGridAdapter(std::shared_ptr<const Embeddings> embeddings)
    : embeddingsPtr_(embeddings)
    , bulkVertexMapper_(embeddings->template gridView<bulkGridId>(), Dune::mcmgVertexLayout())
    {
        // bulk insertion to grid index map
        const auto& bulkGridView = embeddings->template gridView<bulkGridId>();
        bulkInsertionToGridVIdx_.resize(bulkGridView.size(BulkGridView::dimension));
        for (const auto& v : vertices(bulkGridView))
            bulkInsertionToGridVIdx_[embeddings->template insertionIndex<bulkGridId>(v)] = bulkVertexMapper_.index(v);

        // maybe set up index map from hierachy insertion to bulk insertion indices
        makeBulkIndexMap_();

        // determine which bulk vertices lie on facet elements
        bulkVertexIsOnFacetGrid_.resize(bulkGridView.size(BulkGridView::dimension), false);
        const auto& facetGridView = embeddings->template gridView<facetGridId>();
        for (const auto& v : vertices(facetGridView))
        {
            const auto insIdx = embeddings->template insertionIndex<facetGridId>(v);
            const auto highestLevelInsIdx = embeddings->lowDimVertexIndices(facetGridId)[insIdx];
            const auto bulkGridIdx = getBulkGridVertexIndex_(highestLevelInsIdx);
            bulkVertexIsOnFacetGrid_[ bulkGridIdx ] = true;
        }

        // determine the bulk vertex indices that make up facet elements & connectivity
        facetElementCorners_.resize(facetGridView.size(0));
        facetElementsAtBulkVertex_.resize(bulkGridView.size(BulkGridView::dimension));

        std::size_t facetElementCounter = 0;
        for (const auto& element : elements(facetGridView))
        {
            if (isEmbedded(element))
            {
                // obtain the bulk vertex indices of the corners of this element
                const auto numCorners = element.subEntities(FacetGridView::dimension);
                std::vector<BulkIndexType> cornerIndices(numCorners);
                for (int i = 0; i < numCorners; ++i)
                    cornerIndices[i] = bulkGridVertexIndex(element.template subEntity<FacetGridView::dimension>(i));

                // update connectivity map facetVertex -> facetElements
                for (auto bulkVIdx : cornerIndices)
                    facetElementsAtBulkVertex_[bulkVIdx].push_back(facetElementCounter);

                // update facet elements (identified by corners - store them sorted!)
                std::sort(cornerIndices.begin(), cornerIndices.end());
                facetElementCorners_[facetElementCounter] = std::move(cornerIndices);
            }

            facetElementCounter++;
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
        const auto insIdx = embeddingsPtr_->template insertionIndex<facetGridId>(v);
        const auto highestLevelInsIdx = embeddingsPtr_->lowDimVertexIndices(facetGridId)[insIdx];
        return getBulkGridVertexIndex_(highestLevelInsIdx);
    }

    /*!
     * \brief Returns true if the vertex of the d-dimensional grid with
     *        the given vertex index also exists on the (d-1)-dimensional grid
     */
    bool isOnFacetGrid(const BulkGridVertex& v) const
    {
        const auto bulkInsIdx = embeddingsPtr_->template insertionIndex<bulkGridId>(v);
        const auto bulkVIdx = bulkInsertionToGridVIdx_[bulkInsIdx];
        return bulkVertexIsOnFacetGrid_[bulkVIdx];
    }

    /*!
     * \brief Returns true if the given intersection coincides with a facet grid
     *
     * \param element An element of the bulk grid
     * \param intersection An bulk grid intersection within the given element
     */
    bool isOnFacetGrid(const BulkGridElement& element, const BulkGridIntersection& intersection) const
    {
        // Intersection lies on facet grid, if the corners of the intersection make up a facet element
        const auto refElement = BulkReferenceElements::general(element.geometry().type());
        const auto numCorners = intersection.geometry().corners();
        const auto facetIdx = intersection.indexInInside();

        std::vector<BulkIndexType> cornerIndices(numCorners);
        for (int i = 0; i < numCorners; ++i)
            cornerIndices[i] = bulkVertexMapper_.subIndex( element,
                                                           refElement.subEntity(facetIdx, 1, i, BulkGridView::dimension),
                                                           BulkGridView::dimension );

        return composeFacetElement(cornerIndices);
    }

    /*!
     * \brief Returns true if a given set of bulk vertex indices make up a facet grid element
     * \note The vertex indices are NOT the dof indices in the context of models where there
     *       are multiple dofs at one vertex (enriched nodal dofs). Here, we expect the vertex
     *       indices (unique index per vertex).
     */
    template<class IndexStorage>
    bool composeFacetElement(const IndexStorage& bulkVertexIndices) const
    {
        // set up a vector containing all element indices the vertices are connected to
        std::vector<std::size_t> facetElemIndices;
        for (auto bulkVIdx : bulkVertexIndices)
            facetElemIndices.insert( facetElemIndices.end(),
                                     facetElementsAtBulkVertex_[bulkVIdx].begin(),
                                     facetElementsAtBulkVertex_[bulkVIdx].end() );

        // if no facet elements are connected to the vertices this is not on facet grid
        if (facetElemIndices.size() == 0)
            return false;

        // make the container unique
        std::sort(facetElemIndices.begin(), facetElemIndices.end());
        facetElemIndices.erase(std::unique(facetElemIndices.begin(), facetElemIndices.end()), facetElemIndices.end());

        // check if given indices make up a facet element
        auto cornerIndexCopy = bulkVertexIndices;
        std::sort(cornerIndexCopy.begin(), cornerIndexCopy.end());
        for (const auto& facetElemIdx : facetElemIndices)
        {
            const auto& facetElemCorners = facetElementCorners_[facetElemIdx];
            if (facetElemCorners.size() != cornerIndexCopy.size())
                continue;

            if ( std::equal(cornerIndexCopy.begin(), cornerIndexCopy.end(),
                            facetElemCorners.begin(), facetElemCorners.end()) )
                return true;
        }

        // no corresponding facet element found
        return false;
    }

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
    { return embeddingsPtr_->template adjoinedEntityIndices<facetGridId>(e).size(); }

private:
    //! Determine map from the insertion idx of the highest-dimensional grid to bulk insertion idx
    template< bool isHighest = bulkHasHighestDimension, std::enable_if_t<!isHighest, int> = 0 >
    void makeBulkIndexMap_()
    {
        // obtain highest-dimensional grid using the bulkId from grid creator
        const auto& highestLevelGridView = embeddingsPtr_->template gridView<Embeddings::bulkGridId>();
        highestLevelInsertionToBulkInsertionIdx_.resize(highestLevelGridView.size(Embeddings::bulkDim));

        // indices in grid creator stem from the grid file, which has one set of vertices
        // for the entire hierarchy. Thus, the vertex indices we obtain from the grid creator
        // that make up this bulk grid (lower-dimensional within the hierarchy), correspond to
        // the highest-dimensional grid's insertion indices.
        const auto& vertexInsIndices = embeddingsPtr_->lowDimVertexIndices(bulkGridId);
        for (std::size_t insIdx = 0; insIdx < vertexInsIndices.size(); ++insIdx)
            highestLevelInsertionToBulkInsertionIdx_[ vertexInsIndices[insIdx] ] = insIdx;
    }

    //! If the given bulk grid is on the highest level of grid creation, we do not need the map
    template< bool isHighest = bulkHasHighestDimension, std::enable_if_t<isHighest, int> = 0 >
    void makeBulkIndexMap_()
    {}

    //! Obtains the bulk grid vertex index from a given insertion index on the hierarchy
    BulkIndexType getBulkGridVertexIndex_(BulkIndexType highestLevelInsertionIdx) const
    {
        return bulkHasHighestDimension
               ? bulkInsertionToGridVIdx_[ highestLevelInsertionIdx ]
               : bulkInsertionToGridVIdx_[ highestLevelInsertionToBulkInsertionIdx_[ highestLevelInsertionIdx] ];
    }

    // shared pointer to the embedment data
    std::shared_ptr<const Embeddings> embeddingsPtr_;

    // vertex mapper of the bulk grid
    BulkMapper bulkVertexMapper_;

    // data stored on grid vertices
    std::vector<bool> bulkVertexIsOnFacetGrid_;
    std::vector<BulkIndexType> bulkInsertionToGridVIdx_;
    std::vector<BulkIndexType> highestLevelInsertionToBulkInsertionIdx_;

    // data stored for elements on the codim one grid
    std::vector< std::vector<BulkIndexType> > facetElementsAtBulkVertex_;
    std::vector< std::vector<BulkIndexType> > facetElementCorners_;
};

} // end namespace Dumux

#endif
