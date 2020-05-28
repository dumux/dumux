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
 * \copydoc Dumux::EnrichedVertexDofMapper
 */
#ifndef DUMUX_ENRICHED_VERTEX_DOF_MAPPER_HH
#define DUMUX_ENRICHED_VERTEX_DOF_MAPPER_HH

#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include "enrichmenthelper.hh"

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief An indicator class used to mark vertices for enrichment. This
 *        implementation marks all vertices of a given grid of codimension
 *        one for enrichment, except those that are connected to inmersed
 *        boundaries.
 */
class EnrichmentIndicator
{

public:
    /*!
     * \brief Function that marks vertices for enrichment. This implementation
     *        works on the basis of a facet-conforming grid of codimension one
     *        which is used to determine the vertices that should be enriched.
     *        Effectively, all vertices that lie on the given (d-1)-dimensional
     *        grid are marked, except those that are connected to inmersed boundaries.
     *
     * \param vertexMarkers Stores for each vertex if it should be enriched
     * \param gridView The d-dimensional grid for which vertices should be enriched
     * \param vertexMapper Mapper that maps to the vertex indices of the given grid view
     * \param codimOneGridView The view on the (d-1)-dimensional facet-conforming grid
     * \param codimOneGridAdapter Adapter class that allows access to information on the d-
     *                            dimensional grid for entities of the (d-1)-dimensional grid
     */
    template< class GridView,
              class VertexMapper,
              class CodimOneGridView,
              class CodimOneGridAdapter >
    static void markVerticesForEnrichment(std::vector<bool>& vertexMarkers,
                                          const GridView& gridView,
                                          const VertexMapper& vertexMapper,
                                          const CodimOneGridView& codimOneGridView,
                                          const CodimOneGridAdapter& codimOneGridAdapter)
    {
        static constexpr int dim = GridView::dimension;
        static_assert(CodimOneGridView::dimension == dim-1, "Grid dimension mismatch");

        // reset the markers
        vertexMarkers.assign(gridView.size(dim), false);

        // first find all bulk grid vertices on the boundary
        std::vector<bool> isOnBoundary(gridView.size(dim), false);
        for (const auto& e : elements(gridView))
        {
            const auto refElem = referenceElement(e);
            for (const auto& is : intersections(gridView, e))
                if (is.boundary())
                    for (int i = 0; i < is.geometry().corners(); ++i)
                        isOnBoundary[ vertexMapper.subIndex( e,
                                                             refElem.subEntity(is.indexInInside(), 1, i, dim),
                                                             dim ) ] = true;
        }

        // for now, set all markers to true for vertices on codim one grid
        for (const auto& vertex : vertices(codimOneGridView))
            vertexMarkers[codimOneGridAdapter.bulkGridVertexIndex(vertex)] = true;

        // unmark where necessary
        for (const auto& codimOneElement : elements(codimOneGridView))
        {
            // if a codimension one element has less than two embedments  we do not need to enrich
            if (codimOneGridAdapter.numEmbedments(codimOneElement) < 2)
                for (int i = 0; i < codimOneElement.subEntities(dim-1); ++i)
                    vertexMarkers[ codimOneGridAdapter.bulkGridVertexIndex(codimOneElement.template subEntity<dim-1>(i)) ] = false;

            // otherwise, we check for immersed boundaries where we also must not enrich
            else
            {
                const auto refElem = referenceElement(codimOneElement);
                for (const auto& intersection : intersections(codimOneGridView, codimOneElement))
                {
                    // skip if intersection is not on boundary
                    if (!intersection.boundary())
                        continue;

                    // obtain all grid indices of the intersection corners
                    const auto numCorners = intersection.geometry().corners();
                    std::vector<typename GridView::IndexSet::IndexType> vertexIndices(numCorners);
                    for (int i = 0; i < numCorners; ++i)
                    {
                        const auto vIdxLocal = refElem.subEntity(intersection.indexInInside(), 1, i, dim-1);
                        vertexIndices[i] = codimOneGridAdapter.bulkGridVertexIndex( codimOneElement.template subEntity<dim-1>(vIdxLocal) );
                    }

                    // if any of the vertices is on an immersed boudnary, we must not enrich any of them
                    if (std::any_of(vertexIndices.begin(), vertexIndices.end(), [&isOnBoundary] (auto idx) { return !isOnBoundary[idx]; }))
                        std::for_each(vertexIndices.begin(), vertexIndices.end(), [&vertexMarkers] (auto idx) { vertexMarkers[idx] = false; });
                }
            }
        }
    }
};

/*!
 * \ingroup FacetCoupling
 * \brief A vertex mapper that allows for enrichment of nodes. Indication on where to
 *        enrich the nodes is done on the basis of a grid of codimension one living
 *        on the facets of the bulk grid.
 *
 * \tparam GV The Dune::GridView type
 */
template<class GV>
class EnrichedVertexDofMapper
{
    static constexpr int dim = GV::dimension;
    static_assert(dim > 1, "Vertex dof enrichment mapper currently only works for dim > 1!");

    using GIType = typename GV::IndexSet::IndexType;
    using Vertex = typename GV::template Codim<dim>::Entity;
    using Element = typename GV::template Codim<0>::Entity;
    using MCMGMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

public:
    //! export the underlying grid view type
    using GridView = GV;
    //! export the grid index type
    using GridIndexType = GIType;

    //! the constructor
    EnrichedVertexDofMapper(const GV& gridView)
    : gridView_(gridView)
    , elementMapper_(gridView, Dune::mcmgElementLayout())
    , vertexMapper_(gridView, Dune::mcmgVertexLayout())
    {
        initialize_();
    }

    //! constructor taking a layout as additional argument (for compatibility)
    EnrichedVertexDofMapper(const GV& gridView, Dune::MCMGLayout layout)
    : EnrichedVertexDofMapper(gridView)
    {
        if ( !( static_cast<bool>(layout(Dune::GeometryTypes::vertex, dim)) ) )
            DUNE_THROW(Dune::InvalidStateException, "Vertex mapper only makes sense for vertex layout!");
    }

    //! map nodal subentity of codim 0 entity to the grid dof
    GridIndexType subIndex(const Element& e, unsigned int i, unsigned int codim) const
    {
        assert(codim == dim && "Only element corners can be mapped by this mapper");
        return indexMap_[elementMapper_.index(e)][i];
    }

    //! map nodal subentity of codim 0 entity to the grid vertex index
    GridIndexType vertexIndex(const Element& e, unsigned int i, unsigned int codim) const
    {
        assert(codim == dim && "Only element corners can be mapped by this mapper");
        return vertexMapper_.subIndex(e, i, codim);
    }

    //! map nodal entity to the grid vertex index
    GridIndexType vertexIndex(const Vertex& v) const
    {
        assert(Vertex::Geometry::mydimension == 0 && "Only vertices can be mapped by this mapper");
        return vertexMapper_.index(v);
    }

    //! map vertex to the grid dof index
    //! \note This is only valid if there are no enriched vertex dofs!
    //!       We therefore ask this in every call. This means quite some overhead, but
    //!       this mapper is not designed optimally for the case of no enriched nodes.
    template< class EntityType >
    GridIndexType index(const EntityType& e) const
    {
        if (hasEnrichedVertices_)
            DUNE_THROW(Dune::InvalidStateException, "Index map contains enriched vertex dofs. Direct mapping from vertex to index not possible.");

        assert(EntityType::Geometry::mydimension == 0 && "Only vertices can be mapped by this mapper");
        return vertexMapper_.index(e);
    }

    //! returns the number of dofs managed by this mapper
    std::size_t size() const
    { return size_; }

    //! returns true if a vertex dof had been enriched
    bool isEnriched(const Vertex& v)
    { return isEnriched_[ vertexMapper_.index(v) ]; }

    //! the update here simply updates the non-enriched map
    //! enrichment has to be done afterwards!
    void update()
    {
        initialize_();
    }

    /*!
     * \brief Enriches the dof map subject to a (dim-1)-dimensional grid.
     * \note This assumes conforming grids!
     *
     * \param codimOneGridView The grid view of a (dim-1)-dimensional grid conforming
     *                         with the facets of the grid view passed to the constructor,
     *                         indicating on which facets nodal dofs should be enriched.
     * \param codimOneGridAdapter Adapter class that allows access to information on the d-
     *                            dimensional grid for entities of the (d-1)-dimensional grid
     * \param verbose Enable/disable terminal output of the time necessary for enrichment
     */
    template<class CodimOneGridView, class CodimOneGridAdapter>
    void enrich(const CodimOneGridView& codimOneGridView,
                const CodimOneGridAdapter& codimOneGridAdapter,
                bool verbose = false)
    {
        static const int codimOneDim = CodimOneGridView::dimension;
        static_assert(codimOneDim == dim-1, "Grid dimension mismatch!");
        static_assert(codimOneDim == 2 || codimOneDim == 1, "Inadmissible codimension one grid dimension");
        static_assert(int(CodimOneGridView::dimensionworld) == int(GV::dimensionworld), "Grid world dimension mismatch");

        // keep track of time
        Dune::Timer watch;

        // mark vertices for enrichment using the indicator
        EnrichmentIndicator::markVerticesForEnrichment(isEnriched_, gridView_, vertexMapper_, codimOneGridView, codimOneGridAdapter);

        // let the helper class do the enrichment of the index map
        size_ = VertexEnrichmentHelper< GridView, CodimOneGridView >::enrich(indexMap_,
                                                                             isEnriched_,
                                                                             gridView_,
                                                                             vertexMapper_,
                                                                             elementMapper_,
                                                                             codimOneGridView,
                                                                             codimOneGridAdapter);

        // check if new index map contains enriched dofs
        hasEnrichedVertices_ = std::any_of(isEnriched_.begin(), isEnriched_.end(), [] (bool isEnriched) { return isEnriched; });

        if (verbose)
            std::cout << "Vertex dof enrichment took " << watch.elapsed() << " seconds." << std::endl;
    }

private:

    //! initializes the mapper on the basis of the standard Dune mcmgmapper
    void initialize_()
    {
        size_ = gridView_.size(dim);
        hasEnrichedVertices_ = false;
        indexMap_.resize(gridView_.size(0));
        isEnriched_.resize(gridView_.size(dim), false);
        for (const auto& e : elements(gridView_))
        {
            const auto numCorners = e.geometry().corners();
            const auto eIdxGlobal = elementMapper_.index(e);
            indexMap_[eIdxGlobal].resize(numCorners);
            for (unsigned int i = 0; i < numCorners; ++i)
                indexMap_[eIdxGlobal][i] = vertexMapper_.subIndex(e, i, dim);
        }
    }

    // data members
    std::size_t size_;                        //! number of dofs mapped to by this mapper
    const GV gridView_;                       //! the grid view
    MCMGMapper elementMapper_;                //! unmodified element mapper
    MCMGMapper vertexMapper_;                 //! unmodified vertex mapper
    bool hasEnrichedVertices_;                //! keeps track of if vertices are enriched
    std::vector<bool> isEnriched_;            //! keeps track which vertices are enriched
    std::vector< std::vector<GIType> > indexMap_; //! contains the new dof indices
};

} // end namespace Dumux

#endif
