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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \copydoc Dumux::FacetCouplingMapper
 */
#ifndef DUMUX_FACETCOUPLING_MAPPER_HH
#define DUMUX_FACETCOUPLING_MAPPER_HH

#include <dune/common/indices.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Implementation for the coupling mapper that sets up and stores
 *        the coupling maps between two domains of dimension d and (d-1).
 *        The implementations are specific to the discretization method
 *        used in the bulk domain, which is extracted automatically from
 *        the bulk grid geometry. Implementations for the different methods
 *        have to be provided and included at the end of this file.
 *
 * \tparam BulkFVG the d-dimensional finite-volume grid geometry
 * \tparam LowDimFVG the (d-1)-dimensional finite-volume grid geometry
 * \tparam bulkDomainId The domain id of the bulk problem
 * \tparam lowDimDomainId The domain id of the lower-dimensional problem
 * \tparam bulkDM Discretization method used in the bulk domain
 */
template< class BulkFVG,
          class LowDimFVG,
          std::size_t bulkDomainId = 0,
          std::size_t lowDimDomainId = 1,
          DiscretizationMethod bulkDM = BulkFVG::discMethod >
class FacetCouplingMapper;

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Specialization of the mapper class for the case of
 *        three domains with the grid dimensions d, (d-1) & (d-2).
 *
 * \tparam BulkFVG The d-dimensional finite-volume grid geometry
 * \tparam FacetFVG The (d-1)-dimensional finite-volume grid geometry
 * \tparam EdgeFVG The (d-2)-dimensional finite-volume grid geometry
 */
template< class BulkFVG, class FacetFVG, class EdgeFVG,
          std::size_t bulkDomainId = 0,
          std::size_t facetDomainId = 1,
          std::size_t edgeDomainId = 2 >
class FacetCouplingThreeDomainMapper
: public FacetCouplingMapper<BulkFVG, FacetFVG, bulkDomainId, facetDomainId>
, public FacetCouplingMapper<FacetFVG, EdgeFVG, facetDomainId, edgeDomainId>
{
    using BulkFacetMapper = FacetCouplingMapper<BulkFVG, FacetFVG, bulkDomainId, facetDomainId>;
    using FacetEdgeMapper = FacetCouplingMapper<FacetFVG, EdgeFVG, facetDomainId, edgeDomainId>;

public:
    //! Export the coupling stencil type for the provided domain index
    template<std::size_t i>
    using Stencil = typename std::conditional< (i == edgeDomainId),
                                               typename FacetEdgeMapper::template Stencil<i>,
                                               typename BulkFacetMapper::template Stencil<i> >::type;

    //! Export the coupling map type for the provided domain indices
    template<std::size_t i, std::size_t j>
    using CouplingMap = typename std::conditional< (i != edgeDomainId && j != edgeDomainId),
                                                   typename BulkFacetMapper::template CouplingMap<i,j>,
                                                   typename FacetEdgeMapper::template CouplingMap<i,j> >::type;

    /*!
     * \brief Update coupling maps.
     *
     * \param bulkFvGridGeometry The finite-volume grid geometry of the bulk grid
     * \param facetFvGridGeometry The finite-volume grid geometry of the codimension-one grid
     * \param edgeFvGridGeometry The finite-volume grid geometry of the codimension-two grid
     * \param gridCreator Class that contains the grid factories and embedments
     */
    template< class GridCreator >
    void update(const BulkFVG& bulkFvGridGeometry,
                const FacetFVG& facetFvGridGeometry,
                const EdgeFVG& edgeFvGridGeometry,
                const GridCreator& gridCreator)
    {
        BulkFacetMapper::update(bulkFvGridGeometry, facetFvGridGeometry, gridCreator);
        FacetEdgeMapper::update(facetFvGridGeometry, edgeFvGridGeometry, gridCreator);
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

#endif // DUMUX_FACETCOUPLING_COUPLING_MAPPER_HH
