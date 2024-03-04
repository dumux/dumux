// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

/*!
 * \file
 * \ingroup BoxDFMModel
 * \copydoc Dumux::BoxDfmVertexMapper
 */
#ifndef DUMUX_POROUSMEDIUMFLOW_BOXDFM_VERTEX_MAPPER_HH
#define DUMUX_POROUSMEDIUMFLOW_BOXDFM_VERTEX_MAPPER_HH

#include <dumux/multidomain/facet/vertexmapper.hh>

namespace Dumux {

/*!
 * \ingroup BoxDFMModel
 * \brief A vertex mapper that allows for enrichment of dofs on nodes that are located on barrier fractures.
 * \tparam GV The Dune::GridView type
 */
template<class GV>
class BoxDfmVertexMapper : public EnrichedVertexDofMapper<GV>
{
    using ParentType = EnrichedVertexDofMapper<GV>;

    template<class IsBarrier>
    class BarrierIndicator
    {

    public:
        BarrierIndicator(const IsBarrier& isBarrier)
        : isBarrier_{isBarrier}
        , barriersTakePrecendence_(getParam<bool>("BoxDFM.BarriersTakePrecedence", false))
        {}

        /*!
        * \brief marks vertices for enrichment. This implementation
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
        void markVerticesForEnrichment(std::vector<bool>& vertexMarkers,
                                       const GridView& gridView,
                                       const VertexMapper& vertexMapper,
                                       const CodimOneGridView& codimOneGridView,
                                       const CodimOneGridAdapter& codimOneGridAdapter) const
        {
            // first, mark all vertices on the lower-dimensional grid
            EnrichmentIndicator::markVerticesForEnrichment(vertexMarkers, gridView, vertexMapper, codimOneGridView, codimOneGridAdapter);

            static constexpr int dim = GridView::dimension;
            if (barriersTakePrecendence_)
            {
                // find all vertices that are not connected to any barrier elements
                std::vector<bool> touchesBarrier(gridView.size(dim), false);
                for (const auto& codimOneElement : elements(codimOneGridView))
                    if (isBarrier_(codimOneElement))
                        for (int i = 0; i < codimOneElement.subEntities(dim-1); ++i)
                            touchesBarrier[ codimOneGridAdapter.bulkGridVertexIndex(codimOneElement.template subEntity<dim-1>(i)) ] = true;
                for (std::size_t i = 0; i < touchesBarrier.size(); ++i)
                    if (!touchesBarrier[i])
                        vertexMarkers[i] = false;
            }
            else
            {
                // unmark all vertices that are not on barrier elements
                for (const auto& codimOneElement : elements(codimOneGridView))
                    if (!isBarrier_(codimOneElement))
                        for (int i = 0; i < codimOneElement.subEntities(dim-1); ++i)
                            vertexMarkers[ codimOneGridAdapter.bulkGridVertexIndex(codimOneElement.template subEntity<dim-1>(i)) ] = false;
            }
        }

    private:
        const IsBarrier& isBarrier_;
        bool barriersTakePrecendence_;
    };

public:
    using ParentType::ParentType;

    /*!
     * \brief Enriches the dof map subject to a (dim-1)-dimensional grid.
     * \note This assumes conforming grids!
     *
     * \param codimOneGridView The grid view of a (dim-1)-dimensional grid conforming
     *                         with the facets of the grid view passed to the constructor,
     *                         indicating on which facets nodal dofs should be enriched.
     * \param codimOneGridAdapter Adapter class that allows access to information on the d-
     *                            dimensional grid for entities of the (d-1)-dimensional grid
     * \param isBarrier function that takes a (d-1)-dimensional grid element and returns true
                        if it's on a barrier fracture.
     * \param verbose Enable/disable terminal output of the time necessary for enrichment
     */
    template<class CodimOneGridView, class CodimOneGridAdapter, class IsBarrier>
    void enrich(const CodimOneGridView& codimOneGridView,
                const CodimOneGridAdapter& codimOneGridAdapter,
                const IsBarrier& isBarrier,
                bool verbose = false)
    { ParentType::enrich(codimOneGridView, codimOneGridAdapter, BarrierIndicator<IsBarrier>{isBarrier}, verbose); }
};

} // end namespace Dumux

#endif
