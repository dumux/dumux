// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Grid manager specialization for extracting a grid from the facets of a host grid
 */
#ifndef DUMUX_IO_FACET_GRID_MANAGER_HH
#define DUMUX_IO_FACET_GRID_MANAGER_HH

#include <limits>
#include <concepts>
#include <type_traits>
#include <vector>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>

#include "gridmanager.hh"

namespace Dumux {

#ifndef DOXYGEN
namespace FacetGridDetail {

    template<typename Grid, typename HostGridView, typename Selector>
    std::unique_ptr<Grid> makeFacetGrid(const HostGridView& hostGridView, Selector&& selector)
    {
        static constexpr std::size_t undefined = std::numeric_limits<std::size_t>::max();
        static constexpr int domainDim = HostGridView::dimension;

        Dune::MultipleCodimMultipleGeomTypeMapper elementMapper{hostGridView, Dune::mcmgElementLayout()};

        std::vector<unsigned int> localCornerStorage;
        std::vector<std::size_t> domainToFacetVertex(hostGridView.size(domainDim), undefined);

        std::size_t vertexCount = 0;
        Dune::GridFactory<Grid> factory;
        for (const auto& element : elements(hostGridView))
        {
            const auto& refElement = Dune::referenceElement(element);
            const auto& elemGeo = element.geometry();

            for (const auto& is : intersections(hostGridView, element))
            {
                if (!selector(element, is))
                    continue;

                // visit each facet only once
                if (!is.boundary())
                    if (elementMapper.index(is.inside()) > elementMapper.index(is.outside()))
                        continue;

                const auto& isGeo = is.geometry();
                localCornerStorage.clear();
                localCornerStorage.reserve(isGeo.corners());
                for (int c = 0; c < isGeo.corners(); ++c)
                {
                    const auto vIdxLocal = refElement.subEntity(is.indexInInside(), 1, c, domainDim);
                    const auto vIdxGlobal = hostGridView.indexSet().subIndex(element, vIdxLocal, domainDim);
                    if (domainToFacetVertex.at(vIdxGlobal) == undefined)
                    {
                        factory.insertVertex(elemGeo.global(refElement.position(vIdxLocal, domainDim)));
                        domainToFacetVertex[vIdxGlobal] = vertexCount;
                        vertexCount++;
                    }
                    localCornerStorage.push_back(domainToFacetVertex[vIdxGlobal]);
                }

                factory.insertElement(isGeo.type(), localCornerStorage);
            }
        }

        return factory.createGrid();
    }

}  // namespace FacetGridDetail
#endif  // DOXYGEN


namespace Concept {

/*!
 * \ingroup InputOutput
 * \brief Concept for selecting grid intersections to be included in a facet grid.
 */
template<typename T, typename Element, typename Intersection>
concept FacetSelector
    = std::invocable<T, const Element&, const Intersection&>
    and std::convertible_to<bool, std::invoke_result_t<T, const Element&, const Intersection&>>;

} // end namespace Concept

/*!
 * \ingroup InputOutput
 * \brief Grid manager for grids living on the facets of a host grid.
 * \note Each facet is only visited from one side when selecting the facets to be extracted.
 * \note On surface grids as host grids, this implementation assumes that there are no bifurcations.
 */
template<typename HostGrid, typename FacetGrid, typename HostGridManager = GridManager<HostGrid>>
class FacetGridManager
{
    static constexpr int dim = FacetGrid::dimension;
    static constexpr int dimWorld = FacetGrid::dimensionworld;

    static_assert(dim == int(HostGrid::dimension) - 1, "Facet grids must have codimension 1 w.r.t host grid");
    static_assert(dimWorld == int(HostGrid::dimensionworld), "Space dimensions of facet & host grid must match");

    using HostElement = typename HostGrid::template Codim<0>::Entity;
    using HostIntersection = typename HostGrid::LeafGridView::Intersection;

public:
    using Grid = FacetGrid;

    //! Make the grid using an externally created host grid.
    template<Concept::FacetSelector<HostElement, HostIntersection> Selector>
    void init(HostGrid& hostGrid, const Selector& selector)
    {
        facetGrid_ = FacetGridDetail::template makeFacetGrid<Grid>(hostGrid.leafGridView(), selector);
        loadBalance();
    }

    //! Make the grid and create the host grid internally.
    template<Concept::FacetSelector<HostElement, HostIntersection> Selector>
    void init(const Selector& selector, const std::string& paramGroup = "")
    {
        initHostGrid_(paramGroup);
        init(hostGrid_(), selector);
    }

    //! Call loadBalance() function of the grid.
    void loadBalance()
    {
        if (Dune::MPIHelper::getCommunication().size() > 1)
            this->grid().loadBalance();
    }

    //!  Returns a reference to the grid.
    Grid& grid()
    { return *facetGrid_; }

    //1 Returns a const reference to the grid.
    const Grid& grid() const
    { return *facetGrid_; }

    //! Return true if grid data is available
    bool hasGridData() const
    { return false; }

protected:
    void initHostGrid_(const std::string& paramGroup)
    {
        hostGridManager_ = std::make_unique<HostGridManager>();
        hostGridManager_->init(paramGroup);
    }

    HostGrid& hostGrid_()
    { return hostGridManager_->grid(); }

    std::unique_ptr<Grid> facetGrid_{nullptr};
    std::unique_ptr<HostGridManager> hostGridManager_{nullptr};
};


} // namespace Dumux

#endif
