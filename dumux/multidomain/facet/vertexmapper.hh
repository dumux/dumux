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
 * \brief copydoc Dumux::EnrichedVertexDofMapper
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
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief A vertex mapper that allows for enrichment of nodes. Indication on where to
 *        enrich the nodes is done via a lower-dimensional grid living on the facets
 *        of the bulk grid.
 *
 * \tparam GV The Dune::GridView type
 */
template<class GV>
class EnrichedVertexDofMapper
{
    static constexpr int dim = GV::dimension;
    static_assert(dim > 1, "Vertex dof enrichment mapper currently only works for dim > 1!");

    using IT = typename GV::IndexSet::IndexType;
    using Vertex = typename GV::template Codim<dim>::Entity;
    using Element = typename GV::template Codim<0>::Entity;
    using MCMGMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;

    // data members
    std::size_t size_;                        //! number of dofs mapped to by this mapper
    const GV gridView_;                       //! the grid view
    MCMGMapper elementMapper_;                //! unmodified element mapper
    MCMGMapper vertexMapper_;                 //! unmodified vertex mapper
    std::vector< std::vector<IT> > indexMap_; //! contains the new dof indices
    std::vector< bool > isEnriched_;          //! keeps track which vertices are enriched

public:
    //! export the underlying grid view type
    using GridView = GV;
    //! export the grid index type
    using IndexType = IT;

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
    IndexType subIndex(const Element& e, unsigned int i, unsigned int codim) const
    {
        assert(codim == dim && "Only element corners can be mapped by this mapper");
        return indexMap_[elementMapper_.index(e)][i];
    }

    //! map nodal subentity of codim 0 entity to the grid vertex index
    IndexType vertexIndex(const Element& e, unsigned int i, unsigned int codim) const
    {
        assert(codim == dim && "Only element corners can be mapped by this mapper");
        return vertexMapper_.subIndex(e, i, codim);
    }

    //! map nodal subentity of codim 0 entity to the grid vertex index
    IndexType vertexIndex(const Vertex& v) const
    {
        assert(Vertex::Geometry::mydimension == 0 && "Only vertices can be mapped by this mapper");
        return vertexMapper_.index(v);
    }

    //! map vertex to the grid vertex index
    //! \todo TODO this is needed because it is called by the vertex handles.
    //!            Does this work in parallel actually???
    template< class EntityType >
    IndexType index(const EntityType& e) const
    {
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
     * \param facetGridView The grid view of a (dim-1)-dimensional grid conforming
     *                      with the facets of this grid view, indicating on which facets
     *                      nodal dofs should be enriched.
     * \param FacetGridVertexAdapter Allows retrieving the index of (d-1)-dimensional vertices
     *                               within a d-dimensional grid.
     * \param verbose Verbosity level
     */
    template<class FacetGridView, class FacetGridVertexAdapter>
    void enrich(const FacetGridView& facetGridView,
                const FacetGridVertexAdapter& facetGridVertexAdapter,
                bool verbose = false)
    {
        static const int facetDim = FacetGridView::dimension;
        static_assert(facetDim == dim-1, "Grid dimension mismatch!");
        static_assert(facetDim == 2 || facetDim == 1, "Inadmissible facet grid dimension");
        static_assert(int(FacetGridView::dimensionworld) == int(GV::dimensionworld), "Grid world dimension mismatch");

        // keep track of time
        Dune::Timer watch;

        // determine which nodal dofs should actually be enriched
        determineEnrichedVertices_(facetGridView, facetGridVertexAdapter);

        // let the helper class do the enrichment of the index map
        size_ = VertexEnrichmentHelper< GridView, FacetGridView >::enrich(indexMap_,
                                                                          isEnriched_,
                                                                          gridView_,
                                                                          vertexMapper_,
                                                                          elementMapper_,
                                                                          facetGridView,
                                                                          facetGridVertexAdapter);

        if (verbose)
            std::cout << "Vertex dof enrichment took " << watch.elapsed() << " seconds." << std::endl;
    }

private:

    //! initializes the mapper on the basis of the standard Dune mcmgmapper
    void initialize_()
    {
        size_ = gridView_.size(dim);
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

    //! Determines the nodes that should be enriched. This is necessary as dofs on
    //! immersed boundaries of the lower-dimensional grid it should not be enriched.
    template<class FacetGridView, class FacetGridVertexAdapter>
    void determineEnrichedVertices_(const FacetGridView& facetGridView, const FacetGridVertexAdapter& facetGridVertexAdapter)
    {
        // first find all bulk grid vertices on the boundary.
        using ReferenceElements = typename Dune::ReferenceElements<typename GV::ctype, dim>;
        std::vector<bool> isOnDomainBoundary(gridView_.size(dim), false);
        for (const auto& e : elements(gridView_))
        {
            const auto refElem = ReferenceElements::general(e.geometry().type());
            for (const auto& is : intersections(gridView_, e))
                if (is.boundary())
                    for (int i = 0; i < is.geometry().corners(); ++i)
                        isOnDomainBoundary[ vertexMapper_.subIndex( e,
                                                                    refElem.subEntity(is.indexInInside(), 1, i, dim),
                                                                    dim ) ] = true;
        }

        // for now, mark all vertices on facet grid for enrichment
        for (const auto& v : vertices(facetGridView))
            isEnriched_[ facetGridVertexAdapter.bulkGridVertexIndex(v) ] = true;

        // unmark those vertices connected to an immersed boundary
        static const int facetDim = FacetGridView::dimension;
        using FacetReferenceElements = typename Dune::ReferenceElements<typename FacetGridView::ctype, facetDim>;
        for (const auto& facetElement : elements(facetGridView))
        {
            const auto refElem = FacetReferenceElements::general(facetElement.geometry().type());
            for (const auto& facetIs : intersections(facetGridView, facetElement))
            {
                // lambda to obtain the bulk grid vertex indices of the corners
                auto getCornerIndices = [&] (const auto& is)
                {
                    const auto numCorners = is.geometry().corners();
                    std::vector<IndexType> vertexIndices;
                    vertexIndices.reserve(numCorners);
                    for (int i = 0; i < numCorners; ++i)
                    {
                        const auto vIdxLocal = refElem.subEntity(is.indexInInside(), 1, i, facetDim);
                        vertexIndices.push_back( facetGridVertexAdapter.bulkGridVertexIndex(facetElement.template subEntity<facetDim>(vIdxLocal)) );
                    }
                    return vertexIndices;
                };

                // all nodes around non-embedded elements must not be enriched
                if (!facetGridVertexAdapter.isEmbedded(facetElement))
                {
                    const auto vertexIndices = getCornerIndices(facetIs);
                    std::for_each( vertexIndices.begin(),
                                   vertexIndices.end(),
                                   [&] (auto idx) { this->isEnriched_[idx] = false; } );
                }

                // vertices touching immersed boundaries must not be enriched
                else if (facetIs.boundary())
                {
                    const auto vertexIndices = getCornerIndices(facetIs);
                    if ( std::any_of(vertexIndices.begin(),
                                     vertexIndices.end(),
                                     [&isOnDomainBoundary] (auto idx) { return !isOnDomainBoundary[idx]; }) )
                        std::for_each(vertexIndices.begin(),
                                      vertexIndices.end(),
                                      [&] (auto idx) { this->isEnriched_[idx] = false; });
                }
            }
        }
    }
};

} // end namespace Dumux

#endif
