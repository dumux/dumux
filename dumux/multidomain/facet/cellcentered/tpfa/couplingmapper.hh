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
 * \copydoc Dumux::FacetcouplingMapperImplementation
 */
#ifndef DUMUX_FACETCOUPLING_CCTPFA_COUPLING_MAPPER_HH
#define DUMUX_FACETCOUPLING_CCTPFA_COUPLING_MAPPER_HH

#include <dune/common/indices.hh>

#include <dumux/discretization/methods.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmapperbase.hh>

namespace Dumux {

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Base class for the coupling mapper that sets up and stores
 *        the coupling maps between two domains of dimension d and (d-1).
 *        This specialization is for the bulk domain using the cell-centered
 *        scheme with two-point flux approximation.
 *
 * \tparam idOffset Offset added to the mapper-local domain ids for
 *                  the access to the grid quantities in grid creator
 * \tparam BulkFVG The d-dimensional finite-volume grid geometry
 * \tparam LowDimFVG The (d-1)-dimensional finite-volume grid geometry
 */
template<std::size_t idOffset, class BulkFVG, class LowDimFVG>
class FacetCouplingMapperImplementation<idOffset, BulkFVG, LowDimFVG, DiscretizationMethod::cctpfa>
: public virtual FacetCouplingMapperBase<idOffset, BulkFVG, LowDimFVG>
{
    using ParentType = FacetCouplingMapperBase<idOffset, BulkFVG, LowDimFVG>;
    using LowDimElement = typename LowDimFVG::GridView::template Codim<0>::Entity;

    // convenience definitions of domain ids
    static constexpr auto bulkDomainId = Dune::index_constant< idOffset >();
    static constexpr auto lowDimDomainId = Dune::index_constant< idOffset+1 >();

public:
    /*!
     * \brief Update coupling maps. This is the standard
     *        interface required by any mapper implementation.
     *
     * \param bulkFvGridGeometry The finite-volume grid geometry of the bulk grid
     * \param lowDimFvGridGeometry The finite-volume grid geometry of the lower-dimensional grid
     * \param gridCreator Class that contains the grid factories and embedments
     */
    template< class GridCreator >
    void update(const BulkFVG& bulkFvGridGeometry,
                const LowDimFVG& lowDimFvGridGeometry,
                const GridCreator& gridCreator)
    {
        // define the execution policy how to add map entries from an embedment
        auto addEmbedmentPolicy = [&] (auto&& embedments,
                                       const LowDimElement& lowDimElement,
                                       const LowDimFVG& lowDimFvGridGeometry,
                                       const BulkFVG& bulkFvGridGeometry)
        {
            using LowDimIndexType = typename LowDimFVG::GridView::IndexSet::IndexType;
            using BulkIndexType = typename BulkFVG::GridView::IndexSet::IndexType;

            const auto lowDimElemIdx = lowDimFvGridGeometry.elementMapper().index(lowDimElement);
            auto& lowDimData = this->couplingMap_(lowDimDomainId, bulkDomainId)[lowDimElemIdx];

            // find the scvfs in the embedments coinciding with the low dim element
            // since the bulk domain uses tpfa, there is always only going to be one scvf
            for (auto bulkElemIdx : embedments)
            {
                const auto bulkElement = bulkFvGridGeometry.element(bulkElemIdx);

                auto fvGeometry = localView(bulkFvGridGeometry);
                fvGeometry.bindElement(bulkElement);

                bool found = false;
                BulkIndexType embeddedScvfIdx;
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    // for non-boundary faces, it suffices to check if one
                    // of the outside scv indices is in the set of embedments
                    if (!scvf.boundary())
                    {
                        if (std::find(embedments.begin(), embedments.end(), scvf.outsideScvIdx()) != embedments.end())
                        {
                            embeddedScvfIdx = scvf.index();
                            found = true; break;
                        }
                    }

                    // otherwise, do float comparison of element and scvf center
                    else
                    {
                        const auto lowDimGeom = lowDimElement.geometry();
                        const auto eps = lowDimGeom.volume()*1e-8;
                        const auto diffVec = lowDimGeom.center()-scvf.center();

                        using std::abs;
                        if ( std::all_of(diffVec.begin(), diffVec.end(), [eps] (auto coord) { return abs(coord) < eps; }) )
                        {
                            embeddedScvfIdx = scvf.index();
                            found = true; break;
                        }
                    }
                }

                // Error tracking. The boundary scvf detection might has to be improved for very fine grids!?
                if (!found)
                    DUNE_THROW(Dune::InvalidStateException, "Could not find coupling scvf in embedment");

                // add each dof in the low dim element to coupling stencil of the bulk element
                auto& bulkData = this->couplingMap_(bulkDomainId, lowDimDomainId)[bulkElemIdx];
                const auto lowDimElementDofs = LowDimFVG::discMethod == DiscretizationMethod::cctpfa
                                               ? std::vector<LowDimIndexType>( {lowDimElemIdx} )
                                               : this->extractNodalDofs_(lowDimElement, lowDimFvGridGeometry);

                for (auto dofIdx : lowDimElementDofs)
                {
                    bulkData.couplingStencil.push_back( dofIdx );
                    bulkData.couplingScvfs[dofIdx].push_back( embeddedScvfIdx );
                }

                // add embedment (coupling stencil will be done below)
                lowDimData.embedments.emplace_back( bulkElemIdx, std::vector<BulkIndexType>({embeddedScvfIdx}) );
            }

            // embedments = coupling stencil for tpfa
            lowDimData.couplingStencil = std::move(embedments);
        };

        // let the parent do the update subject to the execution policy defined above
        ParentType::update_(bulkFvGridGeometry, lowDimFvGridGeometry, gridCreator, addEmbedmentPolicy);
    }
};

} // end namespace Dumux

#endif // DUMUX_FACETCOUPLING_CCTPFA_COUPLING_MAPPER_HH
