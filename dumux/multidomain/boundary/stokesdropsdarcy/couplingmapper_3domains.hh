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
 * \ingroup DropCoupling
 * \copydoc Dumux::DropCouplingMapper
 */
#ifndef DUMUX_DROPCOUPLING_MAPPER_HH
#define DUMUX_DROPCOUPLING_MAPPER_HH

#include <memory>

#include <dune/common/indices.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup DropCoupling
 * \brief todo
 *
 * \tparam todo
 */
template< class BulkFVG,	// Stokes, Darcy
          class LowDimFVG,	// Interface
          std::size_t bulkId = 0,
          std::size_t lowDimId = 1,
          DiscretizationMethod bulkDM = BulkFVG::discMethod >
class DropCouplingMapper;

/*!
 * \ingroup DropCoupling
 * \brief todo
 *
 * \tparam todo
 */
template< class StokesFVG, class InterfaceFVG, class DarcyFVG,
          std::size_t stokesId = 0,
          std::size_t interfaceId = 1,
          std::size_t darcyId = 2 >
class DropCouplingThreeDomainMapper
: public DropCouplingMapper<StokesFVG, InterfaceFVG, stokesId, interfaceId>
, public DropCouplingMapper<DarcyFVG, InterfaceFVG, darcyId, interfaceId> // TODO order correct?
{
    using StokesDropsMapper = DropCouplingMapper<StokesFVG, InterfaceFVG, stokesId, interfaceId>;
    using DarcyDropsMapper = DropCouplingMapper<DarcyFVG, InterfaceFVG, darcyId, interfaceId>;

    // grid dimensions
    static constexpr int stokesDim = StokesFVG::GridView::dimension;
    static constexpr int interfaceDim = InterfaceFVG::GridView::dimension;
    static constexpr int darcyDim = DarcyFVG::GridView::dimension;

    //! The grid id type
    template<std::size_t id>
    using GridIdType = Dune::index_constant<id>;

public:
    //! export domain ids
    static constexpr auto stokesGridId = Dune::index_constant< stokesId >();
    static constexpr auto interfaceGridId = Dune::index_constant< interfaceId >();
    static constexpr auto darcyGridId = Dune::index_constant< darcyId >();

    // TODO needed for which indices ?!
//    //! Export the coupling stencil type for the provided domain index
//    template<std::size_t i>
//    using Stencil = typename std::conditional< (i == edgeId),
//                                               typename FacetEdgeMapper::template Stencil<i>,
//                                               typename BulkFacetMapper::template Stencil<i> >::type;
//
//    //! Export the coupling map type for the provided domain indices
//    template<std::size_t i, std::size_t j>
//    using CouplingMap = typename std::conditional< (i != edgeId && j != edgeId),
//                                                   typename BulkFacetMapper::template CouplingMap<i,j>,
//                                                   typename FacetEdgeMapper::template CouplingMap<i,j> >::type;

    //! Allow retrievment of grid id for a given grid dimension
    template<int dim>
    static constexpr GridIdType< ( dim == stokesDim ? stokesId : (dim == interfaceDim ? interfaceId : darcyId) ) > gridId()
    { return GridIdType< ( dim == stokesDim ? stokesId : (dim == interfaceDim ? interfaceId : darcyId) ) >(); }

    // TODO necessary??
//    /*!
//     * \brief Update coupling maps.
//     *
//     * \param bulkFvGridGeometry The finite-volume grid geometry of the bulk grid
//     * \param facetFvGridGeometry The finite-volume grid geometry of the codimension-one grid
//     * \param edgeFvGridGeometry The finite-volume grid geometry of the codimension-two grid
//     * \param embeddings Class that contains the embedments among the grids and entity insertion indices
//     */
//    template< class Embeddings >
//    void update(const BulkFVG& bulkFvGridGeometry,
//                const FacetFVG& facetFvGridGeometry,
//                const EdgeFVG& edgeFvGridGeometry,
//                std::shared_ptr<const Embeddings> embeddings)
//    {
//        BulkFacetMapper::update(bulkFvGridGeometry, facetFvGridGeometry, embeddings);
//        FacetEdgeMapper::update(facetFvGridGeometry, edgeFvGridGeometry, embeddings);
//    }

    // TODO necessary??
    //! Pull up the parents' access operators to allow for individual updates
//    using StokesDropsMapper::update;
//    using DarcyDropsMapper::update;
//
//    //! Pull up the parents' access operators
//    using StokesDropsMapper::couplingMap;
//    using DarcyDropsMapper::couplingMap;
};

} // end namespace Dumux

// Here, we have to include all available implementations
// TODO facet or boundary coupling mapper?
#include <dumux/multidomain/facet/box/couplingmapper.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/couplingmapper.hh>

#endif // DUMUX_DROPCOUPLING_COUPLING_MAPPER_HH
