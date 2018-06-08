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
#ifndef DUMUX_FACETCOUPLING_COUPLING_MAPPER_HH
#define DUMUX_FACETCOUPLING_COUPLING_MAPPER_HH

#include <dune/common/indices.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux {

//! Forward declaration of the discretization method-specific implementation between two domains
template<std::size_t idOffset, class BulkFVG, class LowDimFVG, DiscretizationMethod bulkDM>
class FacetCouplingMapperImplementation;

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Class that sets up and stores the coupling maps between an arbitrary
 *        number of sub-domains. These are assumed to be ordered by increasing
 *        codimension. Specializations for the cases of two and three sub-domains
 *        are provided below.
 *
 * \tparam FVG The finite volume grid geometries of the different dimensions
 */
template<class... FVG> class FacetCouplingMapper;

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Class that sets up and stores the coupling maps between two domains
 *        of dimension d and (d-1) in models where the coupling occurs across
 *        the facets of the d-dimensional domain.
 *
 * \note The implemetations are specialized depending on the discretization
 *       method in the bulk domain. We simply inherit from this one here.
 *
 * \tparam BulkFVG The d-dimensional finite-volume grid geometry
 * \tparam LowDimFVG The (d-1)-dimensional finite-volume grid geometry
 */
template<class BulkFVG, class LowDimFVG>
class FacetCouplingMapper<BulkFVG, LowDimFVG>
: public FacetCouplingMapperImplementation<0, BulkFVG, LowDimFVG, BulkFVG::discMethod> {};

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
template<class BulkFVG, class FacetFVG, class EdgeFVG>
class FacetCouplingMapper<BulkFVG, FacetFVG, EdgeFVG>
: public FacetCouplingMapperImplementation<0, BulkFVG, FacetFVG, BulkFVG::discMethod>
, public FacetCouplingMapperImplementation<1, FacetFVG, EdgeFVG, FacetFVG::discMethod>
{
    using BulkFacetMapper = FacetCouplingMapperImplementation<0, BulkFVG, FacetFVG, BulkFVG::discMethod>;
    using FacetEdgeMapper = FacetCouplingMapperImplementation<1, FacetFVG, EdgeFVG, FacetFVG::discMethod>;

public:
    //! Export the coupling stencil type for the provided domain index
    template<std::size_t i>
    using Stencil = typename std::conditional< (i < 2), typename BulkFacetMapper::template Stencil<i>,
                                                        typename FacetEdgeMapper::template Stencil<i> >::type;

    //! Export the coupling map type for the provided domain indices
    template<std::size_t i, std::size_t j>
    using CouplingMap = typename std::conditional< (i < 2 && j < 2), typename BulkFacetMapper::template CouplingMap<i,j>,
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
