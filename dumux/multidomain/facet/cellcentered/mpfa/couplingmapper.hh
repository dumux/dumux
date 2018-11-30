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
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \copydoc Dumux::FacetCouplingMapper
 */
#ifndef DUMUX_CCMPFA_FACETCOUPLING_MAPPER_HH
#define DUMUX_CCMPFA_FACETCOUPLING_MAPPER_HH

#include <dune/common/indices.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmapperbase.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief Base class for the coupling mapper that sets up and stores
 *        the coupling maps between two domains of dimension d and (d-1).
 *        This specialization is for the bulk domain using the cell-centered
 *        scheme with multi-point flux approximation.
 *
 * \tparam BulkFVG The d-dimensional finite-volume grid geometry
 * \tparam LowDimFVG The (d-1)-dimensional finite-volume grid geometry
 * \tparam bulkId The index of the bulk grid within the hierarchy of grids
 * \tparam lowDimId The index of the facet grid within the hierarchy of grids
 */
template<class BulkFVG, class LowDimFVG, std::size_t bulkId, std::size_t lowDimId>
class FacetCouplingMapper<BulkFVG, LowDimFVG, bulkId, lowDimId, DiscretizationMethod::ccmpfa>
: public virtual FacetCouplingMapperBase<BulkFVG, LowDimFVG, bulkId, lowDimId>
{
    using ParentType = FacetCouplingMapperBase<BulkFVG, LowDimFVG, bulkId, lowDimId>;

    using BulkFVElementGeometry = typename BulkFVG::LocalView;
    using BulkSubControlVolumeFace = typename BulkFVG::SubControlVolumeFace;

    using BulkIndexType = typename IndexTraits<typename BulkFVG::GridView>::GridIndex;
    using LowDimIndexType = typename IndexTraits<typename LowDimFVG::GridView>::GridIndex;

    using LowDimElement = typename LowDimFVG::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename LowDimElement::Geometry::GlobalCoordinate;
    static constexpr int lowDimDim = LowDimFVG::GridView::dimension;

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
            const auto lowDimElemIdx = lowDimFvGridGeometry.elementMapper().index(lowDimElement);
            auto& lowDimData = this->couplingMap_(facetGridId, bulkGridId)[lowDimElemIdx];

            // Find the scvfs in the adjoined entities coinciding with the low dim element
            for (auto bulkElemIdx : adjoinedEntityIndices)
            {
                const auto bulkElement = bulkFvGridGeometry.element(bulkElemIdx);

                auto fvGeometry = localView(bulkFvGridGeometry);
                fvGeometry.bind(bulkElement);

                std::vector<BulkIndexType> embeddedScvfIndices;
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    // for non-boundary faces, it suffices to check if one
                    // of the outside scv indices is in the set of embedments
                    if (!scvf.boundary())
                    {
                        if ( std::count(adjoinedEntityIndices.begin(), adjoinedEntityIndices.end(), scvf.outsideScvIdx()) )
                            embeddedScvfIndices.push_back(scvf.index());
                    }

                    // otherwise, do float comparison of element and scvf facet corner
                    else
                    {
                        const auto lowDimGeom = lowDimElement.geometry();
                        const auto eps = lowDimGeom.volume()*1e-8;
                        const auto diffVec = lowDimGeom.center()-scvf.facetCorner();

                        if ( Dune::FloatCmp::eq<GlobalPosition, Dune::FloatCmp::CmpStyle::absolute>(diffVec, GlobalPosition(0.0), eps) )
                            embeddedScvfIndices.push_back(scvf.index());
                    }
                }

                // Error tracking. The boundary scvf detection might has to be improved for very fine grids!?
                if ( embeddedScvfIndices.size() != lowDimElement.subEntities(lowDimDim) )
                    DUNE_THROW(Dune::InvalidStateException, "Could not find all coupling scvfs in embedment");

                // add each dof in the low dim element to coupling stencil of the bulk element and vice versa
                auto& bulkData = this->couplingMap_(bulkGridId, facetGridId)[bulkElemIdx];
                const auto lowDimElementDofs = LowDimFVG::discMethod == DiscretizationMethod::box
                                               ? this->extractNodalDofs_(lowDimElement, lowDimFvGridGeometry)
                                               : std::vector<LowDimIndexType>( {lowDimElemIdx} );

                // add coupling info to all elements/scvfs in interaction volumes
                for (auto dofIdx : lowDimElementDofs)
                    for (auto scvfIdx : embeddedScvfIndices)
                        addCouplingsFromIV_(bulkFvGridGeometry, fvGeometry.scvf(scvfIdx), fvGeometry, lowDimElemIdx, dofIdx);

                // add info on which scvfs coincide with which low dim element
                // TODO order the scvfs according to low dim element corners first !!!
                auto& elemToScvfMap = bulkData.elementToScvfMap[lowDimElemIdx];
                elemToScvfMap.insert(elemToScvfMap.end(), embeddedScvfIndices.begin(), embeddedScvfIndices.end());

                // add embedment
                lowDimData.embedments.emplace_back( bulkElemIdx, std::move(embeddedScvfIndices) );
            }
        };

        // let the parent do the update subject to the execution policy defined above
        ParentType::update_(bulkFvGridGeometry, lowDimFvGridGeometry, embeddings, addCouplingEntryPolicy);

        // lambda to make a container unique
        auto makeUnique = [] (auto& c)
        {
            std::sort(c.begin(), c.end());
            c.erase( std::unique(c.begin(), c.end()), c.end() );
        };

        // lambda to make bulk coupling map entry unique
        auto makeBulkMapEntryUnique = [&makeUnique] (auto& data)
        {
            makeUnique(data.second.couplingStencil);
            makeUnique(data.second.couplingElementStencil);
            std::for_each(data.second.dofToCouplingScvfMap.begin(),
                          data.second.dofToCouplingScvfMap.end(),
                          [&makeUnique] (auto& pair) { makeUnique(pair.second); });
        };

        // make bulk map unique
        auto& bulkCouplingData = this->couplingMap_(bulkGridId, facetGridId);
        std::for_each(bulkCouplingData.begin(), bulkCouplingData.end(), makeBulkMapEntryUnique);

        // coupling stencils might not be unique if box is used in lowdim domain
        if (LowDimFVG::discMethod == DiscretizationMethod::box)
        {
            auto& bulkCouplingData = this->couplingMap_(bulkGridId, facetGridId);
            std::for_each(bulkCouplingData.begin(),
                          bulkCouplingData.end(),
                          [&makeUnique] (auto& pair) { makeUnique(pair.second.couplingStencil); });
        }
    }

private:
    void addCouplingsFromIV_(const BulkFVG& bulkFvGridGeometry,
                             const BulkSubControlVolumeFace& scvf,
                             const BulkFVElementGeometry& fvGeometry,
                             LowDimIndexType lowDimElemIdx,
                             LowDimIndexType lowDimDofIdx)
    {
        const auto& gridIvIndexSets = bulkFvGridGeometry.gridInteractionVolumeIndexSets();
        if (bulkFvGridGeometry.vertexUsesSecondaryInteractionVolume(scvf.vertexIndex()))
            addCouplingsFromIVIndexSet_(gridIvIndexSets.secondaryIndexSet(scvf), fvGeometry, lowDimElemIdx, lowDimDofIdx);
        else
            addCouplingsFromIVIndexSet_(gridIvIndexSets.primaryIndexSet(scvf), fvGeometry, lowDimElemIdx, lowDimDofIdx);
    }

    template< class IVIndexSet >
    void addCouplingsFromIVIndexSet_(const IVIndexSet& indexSet,
                                     const BulkFVElementGeometry& fvGeometry,
                                     LowDimIndexType lowDimElemIdx,
                                     LowDimIndexType lowDimDofIdx)
    {
        auto& lowDimData = this->couplingMap_(facetGridId, bulkGridId)[lowDimElemIdx];
        for (auto bulkElemIdx : indexSet.gridScvIndices())
        {
            auto& bulkData = this->couplingMap_(bulkGridId, facetGridId)[bulkElemIdx];

            lowDimData.couplingStencil.push_back( bulkElemIdx );
            bulkData.couplingStencil.push_back( lowDimDofIdx );
            bulkData.couplingElementStencil.push_back( lowDimElemIdx );
        }

        // insert scvf couplings stemming from this interaction volume
        for (auto scvfIdx : indexSet.gridScvfIndices())
        {
            const auto& scvf = fvGeometry.scvf(scvfIdx);
            auto& bulkData = this->couplingMap_(bulkGridId, facetGridId)[scvf.insideScvIdx()];
            auto& couplingScvfs = bulkData.dofToCouplingScvfMap[lowDimDofIdx];
            couplingScvfs.insert(couplingScvfs.end(), scvfIdx);
        }
    }
};

} // end namespace Dumux

#endif // DUMUX_FACETCOUPLING_CCMPFA_COUPLING_MAPPER_HH
