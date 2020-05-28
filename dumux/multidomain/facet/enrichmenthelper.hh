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
 * \ingroup FacetCoupling
 * \copydoc Dumux::VertexEnrichmentHelper
 */
#ifndef DUMUX_VERTEX_ENRICHMENT_HELPER_HH
#define DUMUX_VERTEX_ENRICHMENT_HELPER_HH

#include <cmath>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/common/math.hh>
#include <dumux/common/indextraits.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Specialization of the enrichment helper class for 2d grids.
 *        In this case, we look for two-dimensional bulk grid elements that
 *        are enclosed by (lie in between) two 1-dimensional facet grid elements.
 *
 * \tparam GridView The grid view of the domain for which nodal dofs should be enriched.
 * \tparam CodimOneGridView The grid view of a (dim-1)-dimensional grid conforming
 *                          with the facets of this grid view, indicating on which facets
 *                          nodal dofs should be enriched.
 */
template< class GridView, class CodimOneGridView >
class VertexEnrichmentHelper
{
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static_assert(dim == 2 || dim == 3, "Grid dimension must be two or three");
    static_assert(dimWorld == int(CodimOneGridView::dimensionworld), "world dimension mismatch");

    using Intersection = typename GridView::Intersection;
    using MCMGMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using ElementPath = std::vector< GridIndexType >;
    using NodalElementPaths = std::vector< std::vector<ElementPath> >;

public:

    /*!
     * \brief Enriches the dof map subject to a (dim-1)-dimensional grid.
     * \note This assumes conforming grids and assumes the index map to be
     *       initialized for the bulk grid already!
     *
     * \param indexMap The index map which is to be updated
     * \param vertexMarkers Markers if vertices should be enriched
     * \param gridView A view on the grid for which vertices should be enriched
     * \param vertexMapper Maps vertices of the given grid view
     * \param elementMapper Maps elements of the given grid view
     * \param codimOneGridView The view on the facet-conforming grid of codim 1
     * \param codimOneGridAdapter Adapter class that allows access to information on the d-
     *                            dimensional grid for entities of the (d-1)-dimensional grid
     * \return the number of dofs after enrichment
     */
    template< class IndexMap, class CodimOneGridAdapter >
    static std::size_t enrich(IndexMap& indexMap,
                              const std::vector<bool>& vertexMarkers,
                              const GridView& gridView,
                              const MCMGMapper& vertexMapper,
                              const MCMGMapper& elementMapper,
                              const CodimOneGridView& codimOneGridView,
                              const CodimOneGridAdapter& codimOneGridAdapter)
    {
        // determine the bulk element paths around vertices marked for enrichment
        NodalElementPaths nodalPaths(gridView.size(dim));
        for (const auto& e : elements(gridView))
        {
            const auto eIdx = elementMapper.index(e);

            std::vector<unsigned int> handledFacets;
            const auto refElement = referenceElement(e);
            for (const auto& is : intersections(gridView, e))
            {
                // do not start path search on boundary intersections
                if (!is.neighbor())
                    continue;

                // if facet has been handled already, skip rest (necesary for e.g. Dune::FoamGrid)
                if (std::find(handledFacets.begin(), handledFacets.end(), is.indexInInside()) != handledFacets.end())
                    continue;

                // first, get indices of all vertices of this intersection
                const auto numCorners = is.geometry().corners();
                std::vector<GridIndexType> faceVertexIndices(numCorners);
                for (int i = 0; i < numCorners; ++i)
                    faceVertexIndices[i] = vertexMapper.subIndex( e,
                                                                  refElement.subEntity(is.indexInInside(), 1, i, dim),
                                                                  dim );

                // start path search on intersections that lie on a facet grid element
                if (codimOneGridAdapter.composeFacetElement(faceVertexIndices))
                {
                    handledFacets.push_back(is.indexInInside());
                    for (int i = 0; i < numCorners; ++i)
                    {
                        // set up path only around those vertices that are to be enriched
                        const auto vIdxGlobal = faceVertexIndices[i];
                        if (!vertexMarkers[vIdxGlobal])
                            continue;

                        // construct the path only if element is not yet contained in another path
                        bool found = false;
                        for (const auto& path : nodalPaths[vIdxGlobal])
                            if (std::find(path.begin(), path.end(), eIdx) != path.end())
                            { found = true; break; }

                        if (!found)
                        {
                            ElementPath path;

                            // Reserve enough memory so that reallocation during the recursive search
                            // does not happen and references are invalidated. Memory is not an issue
                            // here as after enrichment the container is thrown away again.
                            // TODO improve this!?
                            path.reserve(150);
                            path.push_back(eIdx);
                            continuePathSearch_(path, gridView, elementMapper, vertexMapper, codimOneGridAdapter, e, refElement, is, vIdxGlobal);
                            nodalPaths[vIdxGlobal].emplace_back(std::move(path));
                        }
                    }
                }
            }
        }

        // determine the offsets for each bulk vertex index on the basis of the paths found per vertex
        std::vector<std::size_t> bulkVertexIndexOffsets(gridView.size(dim), 0);
        for (const auto& v : vertices(gridView))
        {
            const auto vIdx = vertexMapper.index(v);
            if (vertexMarkers[vIdx])
                bulkVertexIndexOffsets[vIdx] = nodalPaths[vIdx].size()-1;
        }

        // ... and accumulate the offsets
        std::size_t sumOffset = 0;
        std::size_t size = 0;
        for (auto& nodalOffset : bulkVertexIndexOffsets)
        {
            const auto os = nodalOffset;
            nodalOffset = sumOffset;
            sumOffset += os;
            size += (os == 0) ? 1 : os + 1;
        }

        // Now, finally set up the new index map
        for (const auto& e : elements(gridView))
        {
            const auto& eg = e.geometry();
            const auto eIdx = elementMapper.index(e);
            for (int i = 0; i < eg.corners(); ++i)
            {
                const auto origVIdx = vertexMapper.subIndex(e, i, dim);

                // it the node itself is not enriched, simply add offset
                if (!vertexMarkers[origVIdx])
                    indexMap[eIdx][i] += bulkVertexIndexOffsets[origVIdx];

                // find the local index of the path the element belongs to
                else
                {
                    bool found = false;
                    const auto& paths = nodalPaths[ origVIdx ];
                    for (int pathIdx = 0; pathIdx < paths.size(); ++pathIdx)
                    {
                        const auto& curPath = paths[pathIdx];
                        if ( std::find(curPath.begin(), curPath.end(), eIdx) != curPath.end() )
                        {
                            indexMap[eIdx][i] += bulkVertexIndexOffsets[origVIdx] + pathIdx;
                            found = true;
                            break;
                        }
                    }

                    if (!found)
                        DUNE_THROW(Dune::InvalidStateException, "Element not found in any path");
                }
            }
        }

        // return the number of dofs after enrichment
        return size;
    }

private:
    template< class ReferenceElement, class CodimOneGridAdapter >
    static void continuePathSearch_(ElementPath& path,
                                    const GridView& gridView,
                                    const MCMGMapper& elementMapper,
                                    const MCMGMapper& vertexMapper,
                                    const CodimOneGridAdapter& codimOneGridAdapter,
                                    const Element& element,
                                    const ReferenceElement& refElement,
                                    const Intersection& prevIntersection,
                                    GridIndexType vIdxGlobal)
    {
        // on 2d/3d grids, we need to find 1/2 more intersections at vertex
        static constexpr int numIsToFind = dim == 3 ? 2 : 1;

        // keep track of facets handled already while searching
        unsigned int foundCounter = 0;
        std::vector<unsigned int> handledFacets;
        for (const auto& is : intersections(gridView, element))
        {
            // skip if on previous intersection again
            if (is.indexInInside() == prevIntersection.indexInInside())
                continue;

            // skip intersection if handled already (necessary for e.g. Dune::FoamGrid)
            if (std::count(handledFacets.begin(), handledFacets.end(), is.indexInInside()))
                continue;

            // determine all vertex indices of this face
            const auto numCorners = is.geometry().corners();
            std::vector<GridIndexType> faceVertexIndices(numCorners);
            for (int i = 0; i < numCorners; ++i)
                faceVertexIndices[i] = vertexMapper.subIndex( element,
                                                              refElement.subEntity(is.indexInInside(), 1, i, dim),
                                                              dim );

            // we found another intersection at the vertex if it contains the vertex around which we rotate
            if (std::find(faceVertexIndices.begin(), faceVertexIndices.end(), vIdxGlobal) != faceVertexIndices.end())
            {
                foundCounter++;
                handledFacets.push_back(is.indexInInside());

                // keep searching in the outside element only if ...
                // ... this is a not (processor) boundary
                if (!is.neighbor())
                    continue;

                // ... this face does not lie on the facet grid
                if (codimOneGridAdapter.composeFacetElement(faceVertexIndices))
                    continue;

                // ... outside element is not contained yet in path
                const auto outsideElement = is.outside();
                const auto outsideElemIdx = elementMapper.index(outsideElement);
                if (std::find(path.begin(), path.end(), outsideElemIdx) == path.end())
                {
                    const auto idxInOutside = is.indexInOutside();
                    const auto outsideRefElement = referenceElement(outsideElement);
                    path.push_back(outsideElemIdx);

                    // on 2d grids, there is only going to be one more
                    for (const auto& outsideIs : intersections(gridView, outsideElement))
                    {
                        if (outsideIs.indexInInside() == idxInOutside)
                        {
                            // if this is the last intersection to find, return
                            if (foundCounter == numIsToFind)
                                return continuePathSearch_(path,
                                                           gridView,
                                                           elementMapper,
                                                           vertexMapper,
                                                           codimOneGridAdapter,
                                                           outsideElement,
                                                           outsideRefElement,
                                                           outsideIs,
                                                           vIdxGlobal);

                            // otherwise, put one more search on stack
                            else
                                continuePathSearch_(path,
                                                    gridView,
                                                    elementMapper,
                                                    vertexMapper,
                                                    codimOneGridAdapter,
                                                    outsideElement,
                                                    outsideRefElement,
                                                    outsideIs,
                                                    vIdxGlobal);
                        }
                    }
                }
            }
        }

        if (foundCounter != numIsToFind)
            DUNE_THROW(Dune::InvalidStateException, "Found " << foundCounter << " instead of " << numIsToFind << " intersections around vertex");
    }
};

} // end namespace Dumux

#endif
