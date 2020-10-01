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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup FacetCoupling
 * \copydoc Dumux::FacetCouplingMapper
 */
#ifndef DUMUX_FACETCOUPLING_MAPPER_HH
#define DUMUX_FACETCOUPLING_MAPPER_HH

#include <memory>

#include <dune/common/indices.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Implementation for the coupling mapper that sets up and stores
 *        the coupling maps between two domains of dimension d and (d-1).
 *        The implementations are specific to the discretization method
 *        used in the bulk domain, which is extracted automatically from
 *        the bulk grid geometry. Implementations for the different methods
 *        have to be provided and included at the end of this file.
 *
 * \tparam BulkFVG the d-dimensional finite-volume grid geometry
 * \tparam LowDimFVG the (d-1)-dimensional finite-volume grid geometry
 * \tparam bulkId The index of the bulk grid within the hierarchy of grids
 * \tparam lowDimId The index of the facet grid within the hierarchy of grids
 * \tparam bulkDM Discretization method used in the bulk domain
 */
template< class BulkFVG,
          class LowDimFVG,
          std::size_t bulkId = 0,
          std::size_t lowDimId = 1,
          DiscretizationMethod bulkDM = BulkFVG::discMethod >
class FacetCouplingMapper;

/*!
 * \ingroup FacetCoupling
 * \brief Specialization of the mapper class for the case of
 *        three domains with the grid dimensions d, (d-1) & (d-2).
 *
 * \tparam BulkFVG The d-dimensional finite-volume grid geometry
 * \tparam FacetFVG The (d-1)-dimensional finite-volume grid geometry
 * \tparam EdgeFVG The (d-2)-dimensional finite-volume grid geometry
 * \tparam bulkId The index of the bulk grid within the hierarchy of grids
 * \tparam facetId The index of the facet grid within the hierarchy of grids
 * \tparam edgeId The index of the edge grid within the hierarchy of grids
 */
template< class BulkFVG, class FacetFVG, class EdgeFVG,
          std::size_t bulkId = 0,
          std::size_t facetId = 1,
          std::size_t edgeId = 2 >
class FacetCouplingThreeDomainMapper
: public FacetCouplingMapper<BulkFVG, FacetFVG, bulkId, facetId>
, public FacetCouplingMapper<FacetFVG, EdgeFVG, facetId, edgeId>
{
    using BulkFacetMapper = FacetCouplingMapper<BulkFVG, FacetFVG, bulkId, facetId>;
    using FacetEdgeMapper = FacetCouplingMapper<FacetFVG, EdgeFVG, facetId, edgeId>;

    // grid dimensions
    static constexpr int bulkDim = BulkFVG::GridView::dimension;
    static constexpr int facetDim = FacetFVG::GridView::dimension;
    static constexpr int edgeDim = EdgeFVG::GridView::dimension;

    //! The grid id type
    template<std::size_t id>
    using GridIdType = Dune::index_constant<id>;

public:
    //! export domain ids
    static constexpr auto bulkGridId = Dune::index_constant< bulkId >();
    static constexpr auto facetGridId = Dune::index_constant< facetId >();
    static constexpr auto edgeGridId = Dune::index_constant< edgeId >();

    //! Export the coupling stencil type for the provided domain index
    template<std::size_t i>
    using Stencil = typename std::conditional< (i == edgeId),
                                               typename FacetEdgeMapper::template Stencil<i>,
                                               typename BulkFacetMapper::template Stencil<i> >::type;

    //! Export the coupling map type for the provided domain indices
    template<std::size_t i, std::size_t j>
    using CouplingMap = typename std::conditional< (i != edgeId && j != edgeId),
                                                   typename BulkFacetMapper::template CouplingMap<i,j>,
                                                   typename FacetEdgeMapper::template CouplingMap<i,j> >::type;

    //! Allow retrievment of grid id for a given grid dimension
    template<int dim>
    static constexpr GridIdType< ( dim == bulkDim ? bulkId : (dim == facetDim ? facetId : edgeId) ) > gridId()
    { return GridIdType< ( dim == bulkDim ? bulkId : (dim == facetDim ? facetId : edgeId) ) >(); }

    /*!
     * \brief Update coupling maps.
     *
     * \param bulkFvGridGeometry The finite-volume grid geometry of the bulk grid
     * \param facetFvGridGeometry The finite-volume grid geometry of the codimension-one grid
     * \param edgeFvGridGeometry The finite-volume grid geometry of the codimension-two grid
     * \param embeddings Class that contains the embedments among the grids and entity insertion indices
     */
    template< class Embeddings >
    void update(const BulkFVG& bulkFvGridGeometry,
                const FacetFVG& facetFvGridGeometry,
                const EdgeFVG& edgeFvGridGeometry,
                std::shared_ptr<const Embeddings> embeddings)
    {
        BulkFacetMapper::update(bulkFvGridGeometry, facetFvGridGeometry, embeddings);
        FacetEdgeMapper::update(facetFvGridGeometry, edgeFvGridGeometry, embeddings);
    }

    //! Pull up the parents' access operators to allow for individual updates
    using BulkFacetMapper::update;
    using FacetEdgeMapper::update;

    //! Pull up the parents' access operators
    using BulkFacetMapper::couplingMap;
    using FacetEdgeMapper::couplingMap;
};

} // end namespace Dumux

// Here, we have to include all available implementations
#include <dumux/multidomain/facet/box/couplingmapper.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/couplingmapper.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/couplingmapper.hh>

#endif
