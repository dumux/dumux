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
#ifndef DUMUX_CCTPFA_FACETCOUPLING_MAPPER_HH
#define DUMUX_CCTPFA_FACETCOUPLING_MAPPER_HH

#include <dune/common/indices.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmapperbase.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Base class for the coupling mapper that sets up and stores
 *        the coupling maps between two domains of dimension d and (d-1).
 *        This specialization is for the bulk domain using the cell-centered
 *        scheme with two-point flux approximation.
 *
 * \tparam BulkFVG The d-dimensional finite-volume grid geometry
 * \tparam LowDimFVG The (d-1)-dimensional finite-volume grid geometry
 * \tparam bulkId The index of the bulk grid within the hierarchy of grids
 * \tparam lowDimId The index of the facet grid within the hierarchy of grids
 */
template<class BulkFVG, class LowDimFVG, std::size_t bulkId, std::size_t lowDimId>
class FacetCouplingMapper<BulkFVG, LowDimFVG, bulkId, lowDimId, DiscretizationMethod::cctpfa>
: public virtual FacetCouplingMapperBase<BulkFVG, LowDimFVG, bulkId, lowDimId>
{
    using ParentType = FacetCouplingMapperBase<BulkFVG, LowDimFVG, bulkId, lowDimId>;
    using LowDimElement = typename LowDimFVG::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename LowDimElement::Geometry::GlobalCoordinate;

public:
    //! export domain ids
    using ParentType::bulkGridId;
    using ParentType::facetGridId;

    /*!
     * \brief Update coupling maps. This is the standard
     *        interface required by any mapper implementation.
     *
     * \param bulkFvGridGeometry The finite-volume grid geometry of the bulk grid
     * \param lowDimFvGridGeometry The finite-volume grid geometry of the lower-dimensional grid
     * \param embeddings Class that contains the embedments among the grids and entity insertion indices
     */
    template< class Embeddings >
    void update(const BulkFVG& bulkFvGridGeometry,
                const LowDimFVG& lowDimFvGridGeometry,
                std::shared_ptr<const Embeddings> embeddings)
    {
        // define the policy how to add map entries for given lowdim element and adjoined entity indices
        auto addCouplingEntryPolicy = [&] (auto&& adjoinedEntityIndices,
                                           const LowDimElement& lowDimElement,
                                           const LowDimFVG& lowDimFvGridGeometry,
                                           const BulkFVG& bulkFvGridGeometry)
        {
            using LowDimIndexType = typename IndexTraits<typename LowDimFVG::GridView>::GridIndex;
            using BulkIndexType = typename IndexTraits<typename BulkFVG::GridView>::GridIndex;

            const auto lowDimGeometry = lowDimElement.geometry();
            const auto lowDimElemIdx = lowDimFvGridGeometry.elementMapper().index(lowDimElement);
            auto& lowDimData = this->couplingMap_(facetGridId, bulkGridId)[lowDimElemIdx];

            // find the scvfs in the adjoined entities coinciding with the low dim element
            // since the bulk domain uses tpfa, there is always only going to be one scvf
            for (auto bulkElemIdx : adjoinedEntityIndices)
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
                        if ( std::find(adjoinedEntityIndices.begin(),
                                       adjoinedEntityIndices.end(),
                                       scvf.outsideScvIdx()) != adjoinedEntityIndices.end() )
                        {
                            embeddedScvfIdx = scvf.index();
                            found = true; break;
                        }
                    }

                    // otherwise, do float comparison of element and scvf center
                    else
                    {
                        const auto eps = lowDimGeometry.volume()*1e-8;
                        const auto diffVec = lowDimGeometry.center()-scvf.center();

                        if ( Dune::FloatCmp::eq<GlobalPosition, Dune::FloatCmp::CmpStyle::absolute>(diffVec, GlobalPosition(0.0), eps) )
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
                auto& bulkData = this->couplingMap_(bulkGridId, facetGridId)[bulkElemIdx];
                const auto lowDimElementDofs = LowDimFVG::discMethod == DiscretizationMethod::box
                                               ? this->extractNodalDofs_(lowDimElement, lowDimFvGridGeometry)
                                               : std::vector<LowDimIndexType>( {lowDimElemIdx} );

                for (auto dofIdx : lowDimElementDofs)
                {
                    bulkData.couplingStencil.push_back( dofIdx );
                    bulkData.dofToCouplingScvfMap[dofIdx].push_back( embeddedScvfIdx );
                }

                // add info on which scvfs coincide with which low dim element
                bulkData.couplingElementStencil.push_back(lowDimElemIdx);
                bulkData.elementToScvfMap[lowDimElemIdx].push_back( embeddedScvfIdx );

                // add embedment (coupling stencil will be done below)
                lowDimData.embedments.emplace_back( bulkElemIdx, std::vector<BulkIndexType>({embeddedScvfIdx}) );
            }

            // adjoint entity indices = coupling stencil for tpfa
            lowDimData.couplingStencil = std::move(adjoinedEntityIndices);
        };

        // let the parent do the update subject to the execution policy defined above
        ParentType::update_(bulkFvGridGeometry, lowDimFvGridGeometry, embeddings, addCouplingEntryPolicy);

        // coupling stencils might not be unique if box is used in lowdim domain
        if (LowDimFVG::discMethod == DiscretizationMethod::box)
        {
            auto makeStencilUnique = [] (auto& data)
            {
                auto& cs = data.second.couplingStencil;
                std::sort(cs.begin(), cs.end());
                cs.erase( std::unique(cs.begin(), cs.end()), cs.end() );
            };

            auto& bulkCouplingData = this->couplingMap_(bulkGridId, facetGridId);
            std::for_each(bulkCouplingData.begin(), bulkCouplingData.end(), makeStencilUnique);
        }
    }
};

} // end namespace Dumux

#endif
