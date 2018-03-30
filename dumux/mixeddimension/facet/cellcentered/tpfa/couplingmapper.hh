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
 * \copydoc Dumux::CCTpfaFacetCouplingMapper
 */
#ifndef DUMUX_FACETCOUPLING_CCTPFA_COUPLING_MAPPER_HH
#define DUMUX_FACETCOUPLING_CCTPFA_COUPLING_MAPPER_HH

#include <map>
#include <vector>
#include <utility>
#include <algorithm>
#include <type_traits>

#include <dune/common/indices.hh>

namespace Dumux
{
//! Forward declaration of the mapper implementation between two domains
template<std::size_t idOffset, class BulkFVG, class LowDimFVG>
class CCTpfaFacetCouplingTwoDomainMapper;

//! Forward declaration of the mapper class
template<class... FVG> class CCTpfaFacetCouplingMapper;

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Class that sets up and stores the coupling maps between two domains
 *        of dimension d and (d-1) in models where the coupling occurs across
 *        the facets of the d-dimensional domain. This class is to be used for
 *        the cell-centered scheme using tpfa.
 *
 * \tparam BulkFVG the d-dimensional finite-volume grid geometry
 * \tparam LowDimFVG the (d-1)-dimensional finite-volume grid geometry
 */
template<class BulkFVG, class LowDimFVG>
class CCTpfaFacetCouplingMapper<BulkFVG, LowDimFVG>
      : public CCTpfaFacetCouplingTwoDomainMapper<0, BulkFVG, LowDimFVG> {};

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Specialization of the mapper class for the case of
 *        three domains with the grid dimensions d, (d-1) & (d-2).
 *
 * \tparam BulkFVG the d-dimensional finite-volume grid geometry
 * \tparam FacetFVG the (d-1)-dimensional finite-volume grid geometry
 * \tparam EdgeFVG the (d-2)-dimensional finite-volume grid geometry
 */
template<class BulkFVG, class FacetFVG, class EdgeFVG>
class CCTpfaFacetCouplingMapper<BulkFVG, FacetFVG, EdgeFVG>
      : public CCTpfaFacetCouplingTwoDomainMapper<0, BulkFVG, FacetFVG>,
        public CCTpfaFacetCouplingTwoDomainMapper<1, FacetFVG, EdgeFVG>
{
    using BulkFacetMapper = CCTpfaFacetCouplingTwoDomainMapper<0, BulkFVG, FacetFVG>;
    using FacetEdgeMapper = CCTpfaFacetCouplingTwoDomainMapper<1, FacetFVG, EdgeFVG>;

public:
    /*!
     * \brief Update coupling maps
     *
     * \tparam GridCreator Class that contains the grid factories and embedments
     */
    template<class GridCreator>
    void update(const BulkFVG& bulkFvGridGeometry,
                const FacetFVG& facetFvGridGeometry,
                const EdgeFVG& edgeFvGridGeometry,
                const GridCreator& gridCreator)
    {
        BulkFacetMapper::update(bulkFvGridGeometry, facetFvGridGeometry, gridCreator);
        FacetEdgeMapper::update(facetFvGridGeometry, edgeFvGridGeometry, gridCreator);
    }

    using BulkFacetMapper::couplingMap;
    using FacetEdgeMapper::couplingMap;
};

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Implementation that sets up and stores the coupling maps
 *        between two domains of dimension d and (d-1).
 *
 * \tparam idOffset Offset added to the mapper-local domain ids for
 *                  the access to the grid quantities in grid creator
 * \tparam BulkFVG the d-dimensional finite-volume grid geometry
 * \tparam LowDimFVG the (d-1)-dimensional finite-volume grid geometry
 */
template<std::size_t idOffset, class BulkFVG, class LowDimFVG>
class CCTpfaFacetCouplingTwoDomainMapper
{
    using BulkGridView = typename BulkFVG::GridView;
    using BulkElement = typename BulkGridView::template Codim<0>::Entity;
    using BulkIndexType = typename BulkGridView::IndexSet::IndexType;

    using LowDimGridView = typename LowDimFVG::GridView;
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;
    using LowDimIndexType = typename LowDimGridView::IndexSet::IndexType;

    //! make sure the grid geometry combination makes sense
    static constexpr int bulkDim = BulkGridView::dimension;
    static constexpr int lowDimDim = LowDimGridView::dimension;
    static_assert(bulkDim == lowDimDim+1, "Lower-dimensional geometry is not of codim 1 w.r.t. bulk geometry!");

    //! domain ids of the two domains
    static constexpr std::size_t bulkId = idOffset;
    static constexpr std::size_t lowDimId = idOffset+1;

    //! data structure for coupling data bulk->lowDim
    struct BulkCouplingData
    {
        //! indices of the elements coupled to this one
        std::vector<LowDimIndexType> couplingStencil;
        //! for each entry in coupling stencil, we store a list of scvfs whose
        //! fluxes are dependent on that lowDim element (always one face for tpfa)
        std::vector< std::array<BulkIndexType, 1> > couplingScvfs;

        //! returns the list of scvfs whose fluxes are dependent the given element
        const std::array<BulkIndexType, 1>& getCoupledScvfs(LowDimIndexType eIdx) const
        {
            auto it = std::find(couplingStencil.begin(), couplingStencil.end(), eIdx);
            assert(it != couplingStencil.end());
            return couplingScvfs[ std::distance(couplingStencil.begin(), it) ];
        }
    };

    //! data structure for coupling data lowDim->bulk
    struct LowDimCouplingData
    {
        using Embedment = std::pair<BulkIndexType, std::array<BulkIndexType, 1>>;
        //! indices of the elements coupled to this one
        std::vector<BulkIndexType> couplingStencil;
        //! embedments (element + scvf list) of this element
        std::vector<Embedment> embedments;
    };

    //! data structures used for the coupling maps
    using BulkCouplingMap = std::unordered_map<BulkIndexType, BulkCouplingData>;
    using LowDimCouplingMap = std::unordered_map<LowDimIndexType, LowDimCouplingData>;

    //! Couping stencil types to facilitate type export below
    using BulkCouplingStencil = std::vector<LowDimIndexType>;
    using LowDimCouplingStencil = std::vector<BulkIndexType>;

public:
    //! Export the coupling stencil types for the provided domain indices
    template<std::size_t id, std::enable_if_t<(id == bulkId || id == lowDimId), int> = 0>
    using CouplingStencil = typename std::conditional<id == bulkId, BulkCouplingStencil, LowDimCouplingStencil>::type;

    /*!
     * \brief Update coupling maps.
     *
     * \tparam GridCreator Class that contains the grid factories and embedments
     */
    template<class GridCreator>
    void update(const BulkFVG& bulkFvGridGeometry,
                const LowDimFVG& lowDimFvGridGeometry,
                const GridCreator& gridCreator)
    {
        // make sure the provided domain Ids make sense
        using GCBulkGridView = typename GridCreator::template Grid<bulkId>::LeafGridView;
        using GCLowDimGridView = typename GridCreator::template Grid<lowDimId>::LeafGridView;
        static constexpr bool bulkMatch = std::is_same<GCBulkGridView, BulkGridView>::value;
        static constexpr bool lowDimMatch = std::is_same<GCLowDimGridView, LowDimGridView>::value;
        static_assert(bulkMatch, "The bulk domain id does not match the provided bulk grid geometry");
        static_assert(lowDimMatch, "The lowDim domain id does not match the provided lowDim grid geometry");

        // clear & resize
        bulkCouplingData_.clear();
        LowDimCouplingData_.clear();

        // get references on the grid views
        const auto& bulkGridView = bulkFvGridGeometry.gridView();
        const auto& lowDimGridView = lowDimFvGridGeometry.gridView();

        // set up maps between element indices and insertion indices
        std::vector<BulkIndexType> bulkInsertionIdxMap(bulkGridView.size(0));
        for (const auto& e : elements(bulkGridView))
        {
            const auto insIdx = gridCreator.template gridFactory<bulkId>().insertionIndex(e);
            const auto eIdx = bulkFvGridGeometry.elementMapper().index(e);
            bulkInsertionIdxMap[insIdx] = eIdx;
        }

        std::vector<LowDimIndexType> lowDimInsertionIdxMap(lowDimGridView.size(0));
        for (const auto& e : elements(lowDimGridView))
        {
            const auto insIdx = gridCreator.template gridFactory<lowDimId>().insertionIndex(e);
            const auto eIdx = lowDimFvGridGeometry.elementMapper().index(e);
            lowDimInsertionIdxMap[insIdx] = eIdx;
        }

        // set up maps coming from the low dim domain
        for (const auto& element : elements(lowDimFvGridGeometry.gridView()))
        {
            auto embedments = gridCreator.template embedmentEntityIndices<lowDimId>(element);

            // get reference to coupling data for this low dim element
            const auto lowDimIdx = lowDimFvGridGeometry.elementMapper().index(element);
            auto& lowDimData = LowDimCouplingData_[lowDimIdx];

            // turn embedments into actual element indices
            for (unsigned int i = 0; i < embedments.size(); ++i)
                embedments[i] = bulkInsertionIdxMap[ embedments[i] ];

            for (auto bulkIdx : embedments)
            {
                // find the scvf that lies on this low dim element
                const auto bulkElement = bulkFvGridGeometry.element(bulkIdx);

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
                    else
                    {
                        // float comparison of the element and scvf center
                        const auto lowDimGeom = lowDimFvGridGeometry.element(lowDimIdx).geometry();
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

                if (!found) DUNE_THROW(Dune::InvalidStateException, "Could not find coupling scvf in embedment");

                auto& bulkData = bulkCouplingData_[bulkIdx];
                bulkData.couplingStencil.push_back( lowDimIdx );
                bulkData.couplingScvfs.emplace_back( std::array<BulkIndexType,1>({embeddedScvfIdx}) );

                using EB = typename LowDimCouplingData::Embedment;
                lowDimData.couplingStencil.push_back( bulkIdx );
                lowDimData.embedments.emplace_back( EB(bulkIdx, std::array<BulkIndexType,1>({embeddedScvfIdx})) );
            }

            // embedments = coupling stencil for tpfa
            lowDimData.couplingStencil = std::move(embedments);
        }
    }

    //! returns coupling data for bulk->lowDim
    const BulkCouplingMap& couplingMap(Dune::index_constant<bulkId> i,
                                       Dune::index_constant<lowDimId> j) const
    { return bulkCouplingData_; }

    //! returns coupling data for lowDim->bulk
    const LowDimCouplingMap& couplingMap(Dune::index_constant<lowDimId> i,
                                         Dune::index_constant<bulkId> j) const
    { return LowDimCouplingData_; }

private:
    BulkCouplingMap bulkCouplingData_;     //!< stores data on the coupled elements for each bulk element
    LowDimCouplingMap LowDimCouplingData_; //!< stores data on the coupled elements for each low dim element
};

} // end namespace Dumux

#endif // DUMUX_FACETCOUPLING_CCTPFA_COUPLING_MAPPER_HH
