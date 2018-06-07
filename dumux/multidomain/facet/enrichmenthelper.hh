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
 * \brief copydoc Dumux::VertexEnrichmentHelper
 */
#ifndef DUMUX_VERTEX_ENRICHMENT_HELPER_HH
#define DUMUX_VERTEX_ENRICHMENT_HELPER_HH

#include <cmath>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/common/math.hh>

namespace Dumux {

//! Forward declaration of the helper class to perform dof enrichment depending on the grid dimension.
//! Specializations for different facet grid dimensions are provided below.
template< class BulkGridView, class FacetGridView, int facetDim = FacetGridView::dimension >
class VertexEnrichmentHelper;

/*!
 * \brief Specialization of the enrichment helper class for 1d facet grids.
 *        In this case, we look for two-dimensional bulk grid elements that
 *        are enclosed by (lie in between) two 1-dimensional facet grid elements.
 *
 * \tparam BulkGridView The grid view of the bulk domain for which nodal dofs
 *                      should be enriched.
 * \param FacetGridView The grid view of a (dim-1)-dimensional grid conforming
 *                      with the facets of this grid view, indicating on which facets
 *                      nodal dofs should be enriched.
 */
template< class BulkGridView, class FacetGridView >
class VertexEnrichmentHelper< BulkGridView, FacetGridView, 1 >
{
    static constexpr int bulkDim = BulkGridView::dimension;
    static constexpr int facetDim = FacetGridView::dimension;
    static_assert(facetDim == 1, "Grid dimension mismatch");
    static_assert(bulkDim == BulkGridView::dimensionworld, "Enrichment helper still appears to be buggy for surface grids!");

    using BulkIntersection = typename BulkGridView::Intersection;
    using BulkReferenceElements = typename Dune::ReferenceElements<typename BulkGridView::ctype, bulkDim>;
    using BulkMCMGMapper = Dune::MultipleCodimMultipleGeomTypeMapper<BulkGridView>;
    using IndexType = typename BulkGridView::IndexSet::IndexType;
    using BulkElement = typename BulkGridView::template Codim<0>::Entity;
    using GlobalPosition = typename BulkElement::Geometry::GlobalCoordinate;

    using ElementPath = std::vector< IndexType >;
    using NodalElementPaths = std::vector< std::vector<ElementPath> >;

public:

    /*!
     * \brief Enriches the dof map subject to a (dim-1)-dimensional grid.
     * \note This assumes conforming grids and assumes the index map to be
     *       initialized for the bulk grid already!
     *
     * \param facetGridView The grid view of a (dim-1)-dimensional grid conforming
     *                      with the facets of this grid view, indicating on which facets
     *                      nodal dofs should be enriched.
     * \param FacetGridVertexAdapter Allows retrieving the index of (d-1)-dimensional vertices
     *                               within a d-dimensional grid.
     * \return the number of dofs after enrichment
     */
    template< class IndexMap, class FacetGridVertexAdapter >
    static std::size_t enrich(IndexMap& indexMap,
                              const std::vector< bool >& enrichVertex,
                              const BulkGridView& bulkGridView,
                              const BulkMCMGMapper& bulkVertexMapper,
                              const BulkMCMGMapper& bulkElementMapper,
                              const FacetGridView& facetGridView,
                              const FacetGridVertexAdapter& facetGridVertexAdapter)
    {
        // first, find out which bulk vertices lie on the given facet grid
        std::vector< bool > isOnFacetGrid(bulkGridView.size(bulkDim), false);
        for (const auto& v : vertices(facetGridView))
            isOnFacetGrid[ facetGridVertexAdapter.bulkGridVertexIndex(v) ] = true;

        // determine the bulk element paths around vertices marked for enrichment
        NodalElementPaths nodalPaths(bulkGridView.size(bulkDim));
        for (const auto& e : elements(bulkGridView))
        {
            const auto eIdx = bulkElementMapper.index(e);

            std::vector<unsigned int> handledFacets;
            const auto refElement = BulkReferenceElements::general(e.geometry().type());
            for (const auto& is : intersections(bulkGridView, e))
            {
                // do not start path search on boundary intersections
                if (!is.neighbor())
                    continue;

                // if facet has been handled already, skip rest (necesary for e.g. Dune::FoamGrid)
                if (std::find(handledFacets.begin(), handledFacets.end(), is.indexInInside()) != handledFacets.end())
                    continue;

                // first, get indices of all vertices of this intersection
                const auto numCorners = is.geometry().corners();
                std::vector<IndexType> faceVertexIndices(numCorners);
                for (int i = 0; i < numCorners; ++i)
                    faceVertexIndices[i] = bulkVertexMapper.subIndex( e,
                                                                      refElement.subEntity(is.indexInInside(), 1, i, bulkDim),
                                                                      bulkDim );

                // if all vertices of this lie on the facet grid, rotate around the face vertices to find the paths
                if ( std::all_of(faceVertexIndices.begin(),
                                 faceVertexIndices.end(),
                                 [&isOnFacetGrid] (const auto idx) { return isOnFacetGrid[idx]; }) )
                {
                    handledFacets.push_back(is.indexInInside());
                    for (int i = 0; i < numCorners; ++i)
                    {
                        // set up path only around those vertices that are to be enriched
                        const auto vIdxGlobal = faceVertexIndices[i];
                        if (!enrichVertex[vIdxGlobal])
                            continue;

                        // construct the path only if element is not yet contained in one
                        bool found = false;
                        for (const auto& path : nodalPaths[vIdxGlobal])
                            if (std::find(path.begin(), path.end(), eIdx) != path.end())
                            { found = true; break; }

                        if (!found)
                        {
                            ElementPath path;
                            path.push_back(eIdx);
                            continuePathSearch_(path, bulkGridView, bulkElementMapper, bulkVertexMapper, isOnFacetGrid, e, refElement, is, vIdxGlobal);
                            nodalPaths[vIdxGlobal].emplace_back(std::move(path));
                        }
                    }
                }
            }
        }

        // determine the offsets for each bulk vertex index on the basis of the paths found per vertex
        std::vector<std::size_t> bulkVertexIndexOffsets(bulkGridView.size(bulkDim), 0);
        for (const auto& v : vertices(bulkGridView))
        {
            const auto vIdx = bulkVertexMapper.index(v);
            if (enrichVertex[vIdx])
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
        for (const auto& e : elements(bulkGridView))
        {
            const auto& eg = e.geometry();
            const auto eIdx = bulkElementMapper.index(e);
            for (int i = 0; i < eg.corners(); ++i)
            {
                const auto origVIdx = bulkVertexMapper.subIndex(e, i, bulkDim);

                // it the node itself is not enriched, simply add offset
                if (!enrichVertex[origVIdx])
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
    template< class BulkReferenceElement >
    static void continuePathSearch_(ElementPath& path,
                                    const BulkGridView& bulkGridView,
                                    const BulkMCMGMapper& bulkElementMapper,
                                    const BulkMCMGMapper& bulkVertexMapper,
                                    const std::vector<bool>& bulkVertexOnFacetGrid,
                                    const BulkElement& element,
                                    const BulkReferenceElement& refElement,
                                    const BulkIntersection& prevIntersection,
                                    IndexType vIdxGlobal)
    {
        for (const auto& is : intersections(bulkGridView, element))
        {
            // skip the intersection at the vertex found already
            if (is.indexInInside() == prevIntersection.indexInInside())
                continue;

            // determine all vertex indices of this face
            const auto numCorners = is.geometry().corners();
            std::vector<IndexType> faceVertexIndices(numCorners);
            for (int i = 0; i < numCorners; ++i)
                faceVertexIndices[i] = bulkVertexMapper.subIndex( element,
                                                                  refElement.subEntity(is.indexInInside(), 1, i, bulkDim),
                                                                  bulkDim );

            // we found the second intersection if it contains the vertex around which we rotate
            if (std::find(faceVertexIndices.begin(), faceVertexIndices.end(), vIdxGlobal) != faceVertexIndices.end())
            {
                // if this is a (processor) boundary, the search is done
                if (!is.neighbor())
                { return; }

                // if this face lies on the facet grid - we're done too
                if ( std::all_of(faceVertexIndices.begin(),
                                 faceVertexIndices.end(),
                                 [&bulkVertexOnFacetGrid] (const auto idx) { return bulkVertexOnFacetGrid[idx]; }) )
                { return; }

                // otherwise, find common intersection in outside element and proceed
                const auto idxInOutside = is.indexInOutside();
                const auto outsideElement = is.outside();
                const auto outsideRefElement = BulkReferenceElements::general(outsideElement.geometry().type());

                path.push_back( bulkElementMapper.index(outsideElement) );
                for (const auto& outsideIs : intersections(bulkGridView, outsideElement))
                    if (outsideIs.indexInInside() == idxInOutside)
                        return continuePathSearch_(path,
                                                   bulkGridView,
                                                   bulkElementMapper,
                                                   bulkVertexMapper,
                                                   bulkVertexOnFacetGrid,
                                                   outsideElement,
                                                   outsideRefElement,
                                                   outsideIs,
                                                   vIdxGlobal);

                DUNE_THROW(Dune::InvalidStateException, "Could not find intersection in outside element");
            }
        }
    }
};

/*!
 * \brief Specialization of the enrichment helper class for 2d facet grids.
 *        In this case, we find paths of 2d facet grid elements around the facet
 *        grid nodes that enclose different matrix blocks. Afterwards, we check
 *        geometrically which 3d elements are enclosed by which of the paths
 *        and assign the new index accordingly.
 *
 * \tparam BulkGridView The grid view of the bulk domain for which nodal dofs
 *                      should be enriched.
 * \param FacetGridView The grid view of a (dim-1)-dimensional grid conforming
 *                      with the facets of this grid view, indicating on which facets
 *                      nodal dofs should be enriched.
 */
template< class BulkGridView, class FacetGridView >
class VertexEnrichmentHelper< BulkGridView, FacetGridView, 2 >
{
    static constexpr int bulkDim = BulkGridView::dimension;
    static constexpr int facetDim = FacetGridView::dimension;
    static_assert(facetDim == 2, "Grid dimension mismatch");

    using BulkMCMGMapper = Dune::MultipleCodimMultipleGeomTypeMapper<BulkGridView>;
    using FacetMCMGMapper = Dune::MultipleCodimMultipleGeomTypeMapper<FacetGridView>;
    using FacetReferenceElements = typename Dune::ReferenceElements<typename FacetGridView::ctype, facetDim>;

    using IndexType = typename FacetGridView::IndexSet::IndexType;
    using FacetElement = typename FacetGridView::template Codim<0>::Entity;
    using GlobalPosition = typename FacetElement::Geometry::GlobalCoordinate;

    // Structure to store connectivity information around facet grid vertices.
    struct ConnectivityAroundVertex
    {
        // Maps to an adjacent element the two embedded comdim 1 entites touching the vertex
        std::unordered_map< IndexType, Dune::ReservedVector<IndexType, 2> > elementToFacetsMap;
        // stores for each codim 1 entity at the vertex the neighbors connected to it
        std::unordered_map< IndexType, std::unordered_set<IndexType> > facetNeighbors;
    };

    struct ElementPathSegment
    {
        IndexType index;       //!< The index of the cell
        GlobalPosition center; //!< Center of the segment geometry
        GlobalPosition normal; //!< Normal pointing inside the convex hull of element path
    };

    using ElementPath = std::vector< ElementPathSegment >;
    using NodalElementPaths = std::vector< ElementPath >;

public:

    /*!
     * \brief Enriches the dof map subject to a (dim-1)-dimensional grid.
     * \note This assumes conforming grids and assumes the index map to be
     *       initialized for the bulk grid already!
     *
     * \param facetGridView The grid view of a (dim-1)-dimensional grid conforming
     *                      with the facets of this grid view, indicating on which facets
     *                      nodal dofs should be enriched.
     * \param FacetGridVertexAdapter Allows retrieving the index of (d-1)-dimensional vertices
     *                               within a d-dimensional grid.
     *
     * \return the number of dofs after enrichment
     */
    template< class IndexMap, class FacetGridVertexAdapter >
    static std::size_t enrich(IndexMap& indexMap,
                              const std::vector< bool >& enrichVertex,
                              const BulkGridView& bulkGridView,
                              const BulkMCMGMapper& bulkVertexMapper,
                              const BulkMCMGMapper& bulkElementMapper,
                              const FacetGridView& facetGridView,
                              const FacetGridVertexAdapter& facetGridVertexAdapter)
    {
        // determine the element paths around the facet grid vertices
        auto facetVertexMapper = FacetMCMGMapper(facetGridView, Dune::mcmgVertexLayout());
        auto facetElementMapper = FacetMCMGMapper(facetGridView, Dune::mcmgElementLayout());
        auto facetEdgeMapper = FacetMCMGMapper(facetGridView, Dune::mcmgLayout(Dune::Codim<1>()));

        // determine connectivity around vertices
        // Also, for fast access later on, create vectors with element/facet/vertex centers
        std::vector< GlobalPosition > cellCenters(facetGridView.size(0));
        std::vector< GlobalPosition > edgeCenters(facetGridView.size(1));
        std::vector< GlobalPosition > vertexPositions(facetGridView.size(facetDim));
        std::vector< ConnectivityAroundVertex > vertexConnectivities(facetGridView.size(facetDim));

        for (const auto& e : elements(facetGridView))
        {
            // skip non-embedded elements
            if (!facetGridVertexAdapter.isEmbedded(e)) continue;

            const auto eIdx = facetElementMapper.index(e);
            const auto refElem = FacetReferenceElements::general(e.geometry().type());

            cellCenters[eIdx] = e.geometry().center();
            std::vector<unsigned int> handledFacets;
            for (const auto& is : intersections(facetGridView, e))
            {
                const auto idxInInside = is.indexInInside();
                if (std::find(handledFacets.begin(), handledFacets.end(), idxInInside) != handledFacets.end())
                    continue;

                const auto& isGeom = is.geometry();
                const auto edgeIdx = facetEdgeMapper.subIndex(e, idxInInside, 1);
                handledFacets.push_back(idxInInside);
                edgeCenters[edgeIdx] = isGeom.center();
                for (int i = 0; i < isGeom.corners(); ++i)
                {
                    const auto vIdxLocal = refElem.subEntity(idxInInside, 1, i, facetDim);
                    const auto vIdxGlobal = facetVertexMapper.subIndex(e, vIdxLocal, facetDim);
                    const auto& vertex = e.template subEntity<facetDim>(vIdxLocal);
                    const auto idxInBulkGrid = facetGridVertexAdapter.bulkGridVertexIndex(vertex);

                    // Bother setting up the connectivity only for nodes that are enriched
                    if (enrichVertex[idxInBulkGrid])
                    {
                        vertexPositions[vIdxGlobal] = vertex.geometry().center();
                        vertexConnectivities[vIdxGlobal].elementToFacetsMap[eIdx].push_back(edgeIdx);
                        vertexConnectivities[vIdxGlobal].facetNeighbors[edgeIdx].insert(eIdx);
                    }
                }
            }
        }

        // find the paths separating different matrix blocks around the vertex
        const auto numNodes = vertexConnectivities.size();
        std::vector< NodalElementPaths > nodalPaths(numNodes);
        for (const auto& v : vertices(facetGridView))
        {
            if (!enrichVertex[facetGridVertexAdapter.bulkGridVertexIndex(v)])
                continue;

            const auto vIdx = facetVertexMapper.index(v);
            const auto& nodalConnectivity = vertexConnectivities[vIdx];
            std::unordered_set<IndexType> visitedElements;

            for (const auto& elemFacetMap : nodalConnectivity.elementToFacetsMap)
            {
                // only proceed for non-visited elements
                if (visitedElements.count(elemFacetMap.first))
                    continue;

                // set up new path
                ElementPath path;
                const auto firstFacetIdx = elemFacetMap.second[0];
                const auto& vertex = vertexPositions[vIdx];
                const auto& edge = edgeCenters[firstFacetIdx];
                const auto& cell = cellCenters[elemFacetMap.first];
                ElementPathSegment firstSegment;
                firstSegment.index = elemFacetMap.first;
                firstSegment.center = cell;
                firstSegment.normal = crossProduct( edge-cell, vertex-edge );
                firstSegment.normal /= firstSegment.normal.two_norm();

                // recursively add the outside cells to path
                path.push_back(firstSegment);
                bool pathFinished = addOutsideSegmentToPath_(path,
                                                             firstSegment,
                                                             firstFacetIdx,
                                                             cellCenters,
                                                             edgeCenters,
                                                             vertex,
                                                             nodalConnectivity,
                                                             /*invertNormal*/false);

                // maybe continue the path search across the other facet
                if (!pathFinished)
                    addOutsideSegmentToPath_(path,
                                             firstSegment,
                                             elemFacetMap.second[1],
                                             cellCenters,
                                             edgeCenters,
                                             vertex,
                                             nodalConnectivity,
                                             /*invertNormal*/true);

                // make path for the other side of the facet
                ElementPath path2;
                firstSegment.normal *= -1.0;

                // recursively add the outside cells to path
                path2.push_back(firstSegment);
                pathFinished = addOutsideSegmentToPath_(path2,
                                                        firstSegment,
                                                        firstFacetIdx,
                                                        cellCenters,
                                                        edgeCenters,
                                                        vertex,
                                                        nodalConnectivity,
                                                        /*invertNormal*/true);

                // maybe continue path search ...
                if (!pathFinished)
                    addOutsideSegmentToPath_(path2,
                                             firstSegment,
                                             elemFacetMap.second[1],
                                             cellCenters,
                                             edgeCenters,
                                             vertex,
                                             nodalConnectivity,
                                             /*invertNormal*/false);

                // add all cells to set of visited elements
                std::for_each(path.begin(), path.end(), [&] (const auto& segment) { visitedElements.insert(segment.index); });
                std::for_each(path2.begin(), path2.end(), [&] (const auto& segment) { visitedElements.insert(segment.index); });

                // add paths to global container
                nodalPaths[vIdx].emplace_back( std::move(path) );
                nodalPaths[vIdx].emplace_back( std::move(path2) );
            }
        }

        // determine the offsets for each bulk vertex index on the basis of the paths found per vertex
        std::unordered_map<IndexType, typename FacetGridView::IndexSet::IndexType> bulkToFacetIdx;
        std::vector<std::size_t> bulkVertexIndexOffsets(bulkGridView.size(bulkDim), 0);
        for (const auto& v : vertices(facetGridView))
        {
            const auto vIdx = facetVertexMapper.index(v);
            const auto bulkVIdx = facetGridVertexAdapter.bulkGridVertexIndex(v);
            bulkToFacetIdx[bulkVIdx] = vIdx;
            if (enrichVertex[bulkVIdx])
                bulkVertexIndexOffsets[bulkVIdx] = nodalPaths[vIdx].size()-1;
        }

        // ... and accumulate the offsets
        std::size_t size = 0;
        std::size_t sumOffset = 0;
        for (auto& nodalOffset : bulkVertexIndexOffsets)
        {
            const auto os = nodalOffset;
            nodalOffset = sumOffset;
            sumOffset += os;
            size += (os == 0) ? 1 : os + 1;
        }

        // Now, finally set up the new index map
        for (const auto& e : elements(bulkGridView))
        {
            const auto& eg = e.geometry();
            const auto eIdx = bulkElementMapper.index(e);
            for (int i = 0; i < eg.corners(); ++i)
            {
                const auto origVIdx = bulkVertexMapper.subIndex(e, i, bulkDim);

                // it the node itself is not enriched, simply add offset
                if (!enrichVertex[origVIdx])
                    indexMap[eIdx][i] += bulkVertexIndexOffsets[origVIdx];

                // find the path with normals pointing to the element
                else
                {
                    bool found = false;
                    const auto& paths = nodalPaths[ bulkToFacetIdx[origVIdx] ];
                    for (int pathIdx = 0; pathIdx < paths.size(); ++pathIdx)
                    {
                        const auto& curPath = paths[pathIdx];

                        // lambda to test orientation to path segment
                        const auto c = eg.center();
                        auto testOrientation = [&c] (const auto& seg)
                                               { return ( (seg.center-c)*seg.normal ) < 0.0; };

                        // if all normals point towards the element, this is the path
                        if (std::all_of(curPath.begin(), curPath.end(), testOrientation))
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
    //! Adds the next outside segment to an element path around a vertex (facet grid dim = 2)
    static bool addOutsideSegmentToPath_(ElementPath& path,
                                         const ElementPathSegment& lastSegment,
                                         typename FacetGridView::IndexSet::IndexType lastFacetIdx,
                                         const std::vector<GlobalPosition>& cellCenters,
                                         const std::vector<GlobalPosition>& edgeCenters,
                                         const GlobalPosition& vertex,
                                         const ConnectivityAroundVertex& connectivity,
                                         bool invertNormal)
    {
        // list of all neighbors of the current facet
        const auto& neighborList = connectivity.facetNeighbors.find(lastFacetIdx)->second;
        std::vector<typename FacetGridView::IndexSet::IndexType> outsideFacetIndices;
        std::vector<ElementPathSegment> outsideSegments;
        outsideSegments.reserve(neighborList.size()-1);
        outsideFacetIndices.reserve(neighborList.size()-1);

        // make temporary segments for all outside cells
        for (auto eIdx : neighborList)
        {
            if (eIdx != lastSegment.index)
            {
                // if an outside cell is the first of the path, we're done
                if (eIdx == path[0].index)
                    return true;

                // find the other facet (at this node) of this element
                for (auto fIdx : connectivity.elementToFacetsMap.find(eIdx)->second)
                {
                    if (fIdx != lastFacetIdx)
                    {
                        ElementPathSegment seg;
                        seg.index = eIdx;
                        seg.center = cellCenters[eIdx];
                        seg.normal = crossProduct( edgeCenters[fIdx]-seg.center, vertex-edgeCenters[fIdx] );
                        seg.normal /= seg.normal.two_norm();
                        if (invertNormal) seg.normal *= -1.0;
                        outsideSegments.emplace_back( std::move(seg) );
                        outsideFacetIndices.push_back( fIdx );
                        break;
                    }
                }
            }
        }

        // if this is a boundary, return false (rotation not finished)
        if (outsideSegments.size() == 0)
            return false;
        // if it has only one outside segment, the choice is clear
        else if (outsideSegments.size() == 1)
        {
            path.push_back(outsideSegments[0]);
            return addOutsideSegmentToPath_(path,
                                            outsideSegments[0],
                                            outsideFacetIndices[0],
                                            cellCenters,
                                            edgeCenters,
                                            vertex,
                                            connectivity,
                                            invertNormal);
        }
        // if there are more, find the one with the smallest angle
        else
        {
            std::vector<typename GlobalPosition::value_type> angles(outsideSegments.size());

            using std::abs;
            using std::acos;
            static const double pi = acos(-1.0);
            for (int i = 0; i < outsideSegments.size(); ++i)
            {
                const auto& outsideSeg = outsideSegments[i];

                // angle between the normals
                const auto sp = lastSegment.normal*outsideSeg.normal;
                const auto sp2 = lastSegment.normal*(outsideSeg.center-lastSegment.center);
                const auto alpha = acos( sp );

                if ( sp < 1e-7 ) // 0° < alpha < 90° or 270° < alpha < 360°
                {
                    if (sp2 >= 0.0) angles[i] = alpha;
                    else angles[i] = 2*pi-alpha;
                }
                else if ( sp > 1e-7 ) // 90° < alpha < 180° or 180° < alpha < 270°
                {
                    if (sp2 >= 0.0)  angles[i] = pi-alpha;
                    else angles[i] = pi+alpha;
                }
            }

            const auto minAngleIdx = std::distance( angles.begin(), std::min_element(angles.begin(), angles.end()) );
            path.push_back(outsideSegments[minAngleIdx]);
            return addOutsideSegmentToPath_(path,
                                            outsideSegments[minAngleIdx],
                                            outsideFacetIndices[minAngleIdx],
                                            cellCenters,
                                            edgeCenters,
                                            vertex,
                                            connectivity,
                                            invertNormal);
        }
    }
};

} // end namespace Dumux

#endif
