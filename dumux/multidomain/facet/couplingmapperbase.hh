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
 * \copydoc Dumux::FacetCouplingMapperBase
 */
#ifndef DUMUX_FACETCOUPLING_MAPPER_BASE_HH
#define DUMUX_FACETCOUPLING_MAPPER_BASE_HH

#include <vector>
#include <unordered_map>
#include <algorithm>

#include <dune/common/indices.hh>

namespace Dumux {

/*!
 * \ingroup MixedDimension
 * \ingroup MixedDimensionFacet
 * \brief Base class for the coupling mapper that sets up and stores
 *        the coupling maps between two domains of dimension d and (d-1).
 *
 * \tparam BulkFVG the d-dimensional finite-volume grid geometry
 * \tparam LowDimFVG the (d-1)-dimensional finite-volume grid geometry
 * \tparam bulkDomainId The domain id of the bulk problem
 * \tparam lowDimDomainId The domain id of the lower-dimensional problem
 */
template< class BulkFVG,
          class LowDimFVG,
          std::size_t bulkId,
          std::size_t lowDimId>
class FacetCouplingMapperBase
{
    using BulkGridView = typename BulkFVG::GridView;
    using BulkIndexType = typename BulkGridView::IndexSet::IndexType;

    using LowDimGridView = typename LowDimFVG::GridView;
    using LowDimElement = typename LowDimGridView::template Codim<0>::Entity;
    using LowDimIndexType = typename LowDimGridView::IndexSet::IndexType;

    // make sure the grid geometry combination makes sense
    static constexpr int bulkDim = BulkGridView::dimension;
    static constexpr int lowDimDim = LowDimGridView::dimension;
    static_assert(bulkDim == lowDimDim+1, "Lower-dimensional geometry is not of codim 1 w.r.t. bulk geometry!");

    // discretization method of the bulk domain
    static constexpr auto bulkDiscMethod = BulkFVG::discMethod;

    // helper struct to check validity of the given domain id offset
    template< class GridCreator >
    class IsValidGridCreator
    {
        using GCBulkGridView = typename GridCreator::template Grid<bulkId>::LeafGridView;
        using GCLowDimGridView = typename GridCreator::template Grid<lowDimId>::LeafGridView;
        static constexpr bool bulkMatch = std::is_same<GCBulkGridView, BulkGridView>::value;
        static constexpr bool lowDimMatch = std::is_same<GCLowDimGridView, LowDimGridView>::value;
        static_assert(bulkMatch, "The bulk domain id does not match the provided bulk grid geometry");
        static_assert(lowDimMatch, "The lowDim domain id does not match the provided lowDim grid geometry");

    public:
        static constexpr bool value = bulkMatch && lowDimMatch;
    };

    /*!
     * \brief Data structure to store coupling data on the
     *        lower-dimensional grid as seen from the bulk grid.
     */
    struct BulkCouplingData
    {
        //! dof indices of the lower-dimensional domain coupled to this element
        std::vector< LowDimIndexType > couplingStencil;

        //! indices of the elements in the stencil of this bulk element
        std::vector< LowDimIndexType > couplingElementStencil;

        //! For each dof in coupling stencil, we store a list of scvfs whose fluxes depend on it
        std::unordered_map< LowDimIndexType, std::vector<BulkIndexType> > dofToCouplingScvfMap;

        //! For each embedded low dim element, we store the scvfs that coincide with it
        std::unordered_map< LowDimIndexType, std::vector<BulkIndexType> > elementToScvfMap;
    };

    /*!
     * \brief Data structure to store coupling data on the
     *        bulk grid as seen from the lower-dimensional grid.
     */
    struct LowDimCouplingData
    {
        //! indices of the bulk dofs coupled to this element
        std::vector< BulkIndexType > couplingStencil;

        //! for each bulk neighbor we store a list of scvfs that coincide with the low dim element
        using Embedment = std::pair< BulkIndexType, std::vector<BulkIndexType> >;
        std::vector< Embedment > embedments;
    };

    //! data structures used for the coupling maps
    using BulkCouplingMap = std::unordered_map<BulkIndexType, BulkCouplingData>;
    using LowDimCouplingMap = std::unordered_map<LowDimIndexType, LowDimCouplingData>;

    //! Coupling stencil types to facilitate type export below
    using LowDimStencil = std::vector<LowDimIndexType>;
    using BulkStencil = std::vector<BulkIndexType>;

public:
    //! Export the stencil type for the provided domain index
    template<std::size_t id>
    using Stencil = typename std::conditional<id == bulkId, BulkStencil, LowDimStencil>::type;

    //! Export the coupling map type
    template<std::size_t i, std::size_t j>
    using CouplingMap = typename std::conditional<i == bulkId, BulkCouplingMap, LowDimCouplingMap>::type;

    /*!
     * \brief Update coupling maps. This is the standard interface
     *        and has to be overloaded by the implementation.
     */
    template< class GridCreator >
    void update(const BulkFVG& bulkFvGridGeometry,
                const LowDimFVG& lowDimFvGridGeometry,
                const GridCreator& gridCreator)
    { DUNE_THROW(Dune::NotImplemented, "Implementation does not provide an update() function."); }

    //! returns coupling data for bulk -> lowDim
    const BulkCouplingMap& couplingMap(Dune::index_constant<bulkId> bulkDomainId,
                                       Dune::index_constant<lowDimId> lowDimDomainId) const
    { return bulkCouplingData_; }

    //! returns coupling data for lowDim -> bulk
    const LowDimCouplingMap& couplingMap(Dune::index_constant<lowDimId> lowDimDomainId,
                                         Dune::index_constant<bulkId> bulkDomainId) const
    { return lowDimCouplingData_; }

protected:

    /*!
     * \brief Update coupling maps.
     *
     * \param bulkFvGridGeometry The finite-volume grid geometry of the bulk grid
     * \param lowDimFvGridGeometry The finite-volume grid geometry of the lower-dimensional grid
     * \param gridCreator Class that contains the grid factories and embedments
     * \param EmbedmentExecutionPolicy Policy for adding coupling entries for embedments
     */
    template< class GridCreator, typename EmbedmentExecutionPolicy >
    void update_(const BulkFVG& bulkFvGridGeometry,
                 const LowDimFVG& lowDimFvGridGeometry,
                 const GridCreator& gridCreator,
                 EmbedmentExecutionPolicy&& embedmentPolicy)
    {
        // some static assertions on the grid creator
        static_assert(IsValidGridCreator< GridCreator >::value, "Grid type mismatch. Please review the provided domain id offset.");

        // clear data
        bulkCouplingData_.clear();
        lowDimCouplingData_.clear();

        // set up maps between element indices and insertion indices
        const auto& bulkGridFactory = gridCreator.template gridFactory<bulkId>();
        const auto bulkInsertionToElemIdxMap = makeInsertionToGridIndexMap_(bulkGridFactory, bulkFvGridGeometry);

        // set up coupling maps coming from the low dim domain
        for (const auto& element : elements(lowDimFvGridGeometry.gridView()))
        {
            auto embedments = gridCreator.template embedmentEntityIndices<lowDimId>(element);

            // proceed only if embedments were found
            if (embedments.size() == 0)
                continue;

            // turn embedments into actual grid element indices ...
            std::for_each(embedments.begin(), embedments.end(), [&] (auto& idx) { idx = bulkInsertionToElemIdxMap[idx]; });

            // ... and add them
            embedmentPolicy(std::move(embedments),  element, lowDimFvGridGeometry, bulkFvGridGeometry);
        }
    }

    //! Creates a container with the nodal dofs within an element
    template< class FVGridGeometry>
    std::vector< typename FVGridGeometry::GridView::IndexSet::IndexType >
    extractNodalDofs_(const typename FVGridGeometry::GridView::template Codim<0>::Entity& element,
                      const FVGridGeometry& fvGridGeometry)
    {
        static constexpr int dim = FVGridGeometry::GridView::dimension;
        using IndexType = typename FVGridGeometry::GridView::IndexSet::IndexType;

        const auto numCorners = element.subEntities(dim);
        std::vector< IndexType > nodalDofs(numCorners);
        for (int i = 0; i < numCorners; ++i)
            nodalDofs[i] = fvGridGeometry.vertexMapper().subIndex(element, i, dim);

        return nodalDofs;
    }

    //! returns non-const coupling data for bulk -> lowDim
    BulkCouplingMap& couplingMap_(Dune::index_constant<bulkId> bulkDomainId,
                                  Dune::index_constant<lowDimId> lowDimDomainId)
    { return bulkCouplingData_; }

    //! returns non-const coupling data for lowDim -> bulk
    LowDimCouplingMap& couplingMap_(Dune::index_constant<lowDimId> lowDimDomainId,
                                    Dune::index_constant<bulkId> bulkDomainId)
    { return lowDimCouplingData_; }

private:

    //! Creates the map from element insertion index to grid element index
    template< class GridFactory, class FVGridGeometry>
    std::vector< typename FVGridGeometry::GridView::IndexSet::IndexType >
    makeInsertionToGridIndexMap_(const GridFactory& gridFactory, const FVGridGeometry& fvGridGeometry) const
    {
        using IndexType = typename FVGridGeometry::GridView::IndexSet::IndexType;

        std::vector< IndexType > map(fvGridGeometry.gridView().size(0));
        for (const auto& e : elements(fvGridGeometry.gridView()))
            map[ gridFactory.insertionIndex(e) ] = fvGridGeometry.elementMapper().index(e);

        return map;
    }

    BulkCouplingMap bulkCouplingData_;     //!< stores data on the coupled elements for each bulk element
    LowDimCouplingMap lowDimCouplingData_; //!< stores data on the coupled elements for each low dim element
};

} // end namespace Dumux

#endif // DUMUX_FACETCOUPLING_COUPLING_MAPPER_BASE_HH
