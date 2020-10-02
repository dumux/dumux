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
#ifndef DUMUX_BOX_FACETCOUPLING_MAPPER_HH
#define DUMUX_BOX_FACETCOUPLING_MAPPER_HH

#include <dune/common/indices.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmapperbase.hh>
#include <dumux/multidomain/facet/codimonegridadapter.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Base class for the coupling mapper that sets up and stores
 *        the coupling maps between two domains of dimension d and (d-1).
 *        This specialization is for the bulk domain using the box scheme.
 *
 * \tparam BulkFVG The d-dimensional finite-volume grid geometry
 * \tparam LowDimFVG The (d-1)-dimensional finite-volume grid geometry
 * \tparam bulkId The domain id of the bulk problem
 * \tparam lowDimId The domain id of the lower-dimensional problem
 */
template<class BulkFVG, class LowDimFVG, std::size_t bulkId, std::size_t lowDimId>
class FacetCouplingMapper<BulkFVG, LowDimFVG, bulkId, lowDimId, DiscretizationMethod::box>
: public virtual FacetCouplingMapperBase<BulkFVG, LowDimFVG, bulkId, lowDimId>
{
    using ParentType = FacetCouplingMapperBase<BulkFVG, LowDimFVG, bulkId, lowDimId>;
    using LowDimElement = typename LowDimFVG::GridView::template Codim<0>::Entity;

    // dimensions of the two grids
    using BulkGridView = typename BulkFVG::GridView;
    using LowDimGridView = typename LowDimFVG::GridView;
    static constexpr int bulkDim = BulkGridView::dimension;
    static constexpr int lowDimDim = LowDimGridView::dimension;

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
        // forward to update function with instantiated vertex adapter
        using GridAdapter = CodimOneGridAdapter<Embeddings, bulkGridId, facetGridId>;
        update(bulkFvGridGeometry, lowDimFvGridGeometry, embeddings, GridAdapter(embeddings));
    }

    /*!
     * \brief Update coupling maps with a given grid adapter.
     *
     * \param bulkFvGridGeometry The finite-volume grid geometry of the bulk grid
     * \param lowDimFvGridGeometry The finite-volume grid geometry of the lower-dimensional grid
     * \param embeddings Class that contains the embedments among the grids and entity insertion indices
     * \param codimOneGridAdapter Allows direct access to data on the bulk grid for lowdim grid entities
     */
    template< class Embeddings, class CodimOneGridAdapter >
    void update(const BulkFVG& bulkFvGridGeometry,
                const LowDimFVG& lowDimFvGridGeometry,
                std::shared_ptr<const Embeddings> embeddings,
                const CodimOneGridAdapter& codimOneGridAdapter)
    {
        // define the policy how to add map entries for given lowdim element and adjoined entity indices
        auto addCouplingEntryPolicy = [&] (auto&& adjoinedEntityIndices,
                                           const LowDimElement& lowDimElement,
                                           const LowDimFVG& lowDimFvGridGeometry,
                                           const BulkFVG& bulkFvGridGeometry)
        {
            using LowDimIndexType = typename IndexTraits<LowDimGridView>::GridIndex;
            using BulkIndexType = typename IndexTraits<BulkGridView>::GridIndex;

            const auto lowDimElemIdx = lowDimFvGridGeometry.elementMapper().index(lowDimElement);
            auto& lowDimData = this->couplingMap_(facetGridId, bulkGridId)[lowDimElemIdx];

            // determine corner indices (in bulk grid indices)
            const auto& eg = lowDimElement.geometry();
            const auto numElementCorners = eg.corners();
            std::vector<BulkIndexType> elementCorners(numElementCorners);
            for (int i = 0; i < numElementCorners; ++i)
                elementCorners[i] = codimOneGridAdapter.bulkGridVertexIndex(lowDimElement.template subEntity<lowDimDim>(i));

            // save unsorted set of corner indices and search scvfs in adjoined entities
            const auto unsortedElemCorners = elementCorners;
            std::sort(elementCorners.begin(), elementCorners.end());
            for (auto bulkElemIdx : adjoinedEntityIndices)
            {
                const auto bulkElement = bulkFvGridGeometry.element(bulkElemIdx);
                const auto bulkRefElem = referenceElement(bulkElement);

                // find the bulk element facet that lies on this low dim element (assumes conformity!)
                bool found = false;
                unsigned int coupledFacetIndex;
                std::vector<unsigned int> handledFacets;
                for (const auto& is : intersections(bulkFvGridGeometry.gridView(), bulkElement))
                {
                    // skip already handled facets (necessary for e.g. Dune::FoamGrid)
                    if (std::count(handledFacets.begin(), handledFacets.end(), is.indexInInside()))
                        continue;

                    handledFacets.push_back(is.indexInInside());

                    // determine if it lies on low dim element by comparing corner indices
                    const auto numCorners = is.geometry().corners();
                    std::vector<BulkIndexType> facetIndices(numCorners);
                    for (int i = 0; i < numCorners; ++i)
                    {
                        const auto vIdxLocal = bulkRefElem.subEntity(is.indexInInside(), 1, i, bulkDim);
                        facetIndices[i] = bulkFvGridGeometry.vertexMapper().vertexIndex(bulkElement.template subEntity<bulkDim>(vIdxLocal));
                    }

                    std::sort(facetIndices.begin(), facetIndices.end());
                    if ( std::equal(facetIndices.begin(), facetIndices.end(), elementCorners.begin(), elementCorners.end()) )
                    {
                        coupledFacetIndex = is.indexInInside();
                        found = true; break;
                    }
                }

                // ensure that we found the facet!
                if (!found)
                    DUNE_THROW(Dune::InvalidStateException, "Could not find the bulk element coupling facet!");

                // we should always find numElementCorners coupling scvfs
                auto fvGeometry = localView(bulkFvGridGeometry);
                fvGeometry.bindElement(bulkElement);

                unsigned int foundCounter = 0;
                std::vector<BulkIndexType> embeddedScvfIndices(numElementCorners);
                for (const auto& scvf : scvfs(fvGeometry))
                {
                    if (scvf.interiorBoundary() && scvf.facetIndexInElement() == coupledFacetIndex)
                    {
                        // we want to order the set of scvfs lying on the lower-dimensional element such that the i-th scvf
                        // coincides with the i-th low dim element corner. This will help later to identify which scvf fluxes
                        // enter which scv of the low dim element if the lower-dimensional domain uses the box scheme
                        const auto vIdxLocal = bulkRefElem.subEntity(coupledFacetIndex, 1, scvf.indexInElementFacet(), bulkDim);
                        const auto vIdxGlobal = bulkFvGridGeometry.vertexMapper().vertexIndex(bulkElement, vIdxLocal, bulkDim);
                        const auto it = std::find(unsortedElemCorners.begin(), unsortedElemCorners.end(), vIdxGlobal);
                        assert(it != unsortedElemCorners.end());
                        const auto lowDimElemLocalCornerIdx = std::distance(unsortedElemCorners.begin(), it);
                        embeddedScvfIndices[lowDimElemLocalCornerIdx] = scvf.index();
                        foundCounter++;
                    }
                }

                // ensure we found all scvfs
                if (foundCounter != numElementCorners)
                    DUNE_THROW(Dune::InvalidStateException, "Found " << foundCounter << " instead of " << numElementCorners << " coupling scvfs in the bulk element");

                // add each dof in the low dim element to coupling stencil of the bulk element
                auto& bulkData = this->couplingMap_(bulkGridId, facetGridId)[bulkElemIdx];
                const auto lowDimElementDofs = LowDimFVG::discMethod == DiscretizationMethod::box
                                               ? this->extractNodalDofs_(lowDimElement, lowDimFvGridGeometry)
                                               : std::vector<LowDimIndexType>({lowDimElemIdx});

                for (auto lowDimDofIdx : lowDimElementDofs)
                {
                    bulkData.couplingStencil.push_back(lowDimDofIdx);
                    auto& couplingScvfs = bulkData.dofToCouplingScvfMap[lowDimDofIdx];
                    couplingScvfs.insert(couplingScvfs.end(), embeddedScvfIndices.begin(), embeddedScvfIndices.end());
                }

                // add info on which scvfs coincide with which low dim element
                bulkData.couplingElementStencil.push_back(lowDimElemIdx);
                auto& elemToScvfMap = bulkData.elementToScvfMap[lowDimElemIdx];
                elemToScvfMap.insert(elemToScvfMap.end(), embeddedScvfIndices.begin(), embeddedScvfIndices.end());

                // add embedment and the bulk cell dofs to coupling stencil of low dim element
                lowDimData.embedments.emplace_back(bulkElemIdx, std::move(embeddedScvfIndices));

                const auto bulkElementDofs = this->extractNodalDofs_(bulkElement, bulkFvGridGeometry);
                for (auto bulkDofIdx : bulkElementDofs)
                    lowDimData.couplingStencil.push_back(bulkDofIdx);
            }
        };

        // let the parent do the update subject to the execution policy defined above
        ParentType::update_(bulkFvGridGeometry, lowDimFvGridGeometry, embeddings, addCouplingEntryPolicy);

        // coupling stencils might not be unique with the policy above
        auto makeStencilUnique = [] (auto& data)
        {
            auto& cs = data.second.couplingStencil;
            std::sort(cs.begin(), cs.end());
            cs.erase( std::unique(cs.begin(), cs.end()), cs.end() );
        };

        auto& lowDimCouplingData = this->couplingMap_(facetGridId, bulkGridId);
        std::for_each(lowDimCouplingData.begin(), lowDimCouplingData.end(), makeStencilUnique);

        // bulk coupling stencil is only non-unique if box is used
        if (LowDimFVG::discMethod == DiscretizationMethod::box)
        {
            auto& bulkCouplingData = this->couplingMap_(bulkGridId, facetGridId);
            std::for_each(bulkCouplingData.begin(), bulkCouplingData.end(), makeStencilUnique);
        }
    }
};

} // end namespace Dumux

#endif
