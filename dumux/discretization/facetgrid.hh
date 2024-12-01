// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Class to extract facets of a grid into a surface grid.
 */
#ifndef DUMUX_DISCRETIZATION_FACET_GRID_HH
#define DUMUX_DISCRETIZATION_FACET_GRID_HH

#include <config.h>

#include <vector>
#include <memory>
#include <utility>
#include <concepts>
#include <optional>
#include <algorithm>
#include <type_traits>
#include <unordered_map>
#include <limits>
#include <ranges>

#include <dune/common/reservedvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dumux/geometry/geometryintersection.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace FacetGridDetail {

    template<typename Geo1, typename Geo2>
    bool intersect(const Geo1& geo1, const Geo2& geo2)
    {
        using Algo = GeometryIntersection<Geo1, Geo2>;
        typename Algo::Intersection r;
        return Algo::intersection(geo1, geo2, r);
    }

    template<typename GridView>
    auto makeBoundingBoxTree(const GridView& gridView)
    {
        using EntitySet = GridViewGeometricEntitySet<GridView>;
        return BoundingBoxTree<EntitySet>{std::make_shared<EntitySet>(gridView)};
    }

    template<typename ScvfGeo, typename FacetGridView, typename FacetBoundingBoxTree>
    std::optional<std::size_t> coupledElementIndex(const ScvfGeo& scvfGeo,
                                                   const FacetGridView& facetGridView,
                                                   const FacetBoundingBoxTree& bboxTree)
    {
        // TODO: use `intersectingEntities` once it has been made robust for bboxes with width=0 in one direction
        for (const auto& element : elements(facetGridView))
            if (intersect(scvfGeo, element.geometry()))
                return bboxTree.entitySet().index(element);
        return {};
    }

    template<typename Grid, typename GridGeometry, typename  Selector>
    std::pair<std::unique_ptr<Grid>, std::vector<std::size_t>>
    makeFacetGrid(const GridGeometry& gridGeometry, Selector&& selector)
    {
        static constexpr std::size_t undefined = std::numeric_limits<std::size_t>::max();
        static constexpr int domainDim = GridGeometry::GridView::dimension;
        static_assert(int(Grid::dimension) == domainDim - 1);

        const auto& gridView = gridGeometry.gridView();
        const auto& vertexMapper = gridGeometry.vertexMapper();
        const auto& elementMapper = gridGeometry.elementMapper();

        std::vector<unsigned int> localCornerStorage;
        std::vector<std::size_t> domainToFacetVertex(gridView.size(domainDim), undefined);

        std::size_t vertexCount = 0;
        Dune::GridFactory<Grid> factory;
        for (const auto& element : elements(gridView))
        {
            const auto& refElement = Dune::referenceElement(element);
            const auto& elemGeo = element.geometry();

            for (const auto& is : intersections(gridView, element))
            {
                if (!selector(is))
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
                    const auto vIdxGlobal = vertexMapper.subIndex(element, vIdxLocal, domainDim);
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

        auto grid = factory.createGrid();

        // update map from insertion to grid vertex indices
        std::vector<std::size_t> insertionToGridIndex(grid->leafGridView().size(Grid::dimension), 0);
        for (const auto& vertex : vertices(grid->leafGridView()))
            insertionToGridIndex[factory.insertionIndex(vertex)] = grid->leafGridView().indexSet().index(vertex);

        std::size_t domainVertexIndex = 0;
        std::vector<std::size_t> facetToDomainVertexMap(grid->leafGridView().size(Grid::dimension));
        for (auto& facetVertexIndex : domainToFacetVertex)
        {
            if (facetVertexIndex != undefined)
                facetToDomainVertexMap[insertionToGridIndex[facetVertexIndex]] = domainVertexIndex;
            domainVertexIndex++;
        }

        return std::make_pair(std::move(grid), std::move(facetToDomainVertexMap));
    }

}  // namespace FacetGridDetail
#endif  // DOXYGEN

template<typename T, typename GridGeometry>
concept FacetSelectorFor = std::invocable<T, const typename GridGeometry::GridView::Intersection&>
    and std::is_convertible_v<std::invoke_result_t<T, const typename GridGeometry::GridView::Intersection&>, bool>;

/*!
 * \ingroup Discretization
 * \brief Data structure to hold a grid that consists of facets of a discretization.
 * \tparam Grid The facet grid type
 * \tparam GG The discretization from which the facet grid has been extracted.
 * \note If the grid geometry is a surface grid with bifurcations, the result depends on the grid implementation.
 *       On surface grids, this class can only be safely used if there are no bifurcations between the selected facets.
 */
template<typename Grid, typename GG>
class FacetGrid
{
 public:
    using DomainGridGeometry = GG;
    using GridView = typename Grid::LeafGridView;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

    //! Construct a facet grid from a grid geometry and a facet selector
    template<FacetSelectorFor<DomainGridGeometry> Selector>
    explicit FacetGrid(std::shared_ptr<const DomainGridGeometry> gridGeometry, Selector&& selector)
    : FacetGrid(gridGeometry, FacetGridDetail::makeFacetGrid<Grid>(*gridGeometry, std::forward<Selector>(selector)))
    {}

    //! Return a view on the facet grid
    GridView gridView() const
    { return facetGrid_->leafGridView(); }

    //! Return the domain from which this facet grid was extracted
    std::shared_ptr<const DomainGridGeometry> domainGridGeometry() const
    { return domain_; }

    //! Return the index of the given vertex within the domain
    std::size_t domainVertexIndexOf(const Vertex& v) const
    { return facetToDomainVertexIndex_.at(gridView().indexSet().index(v)); }

 private:
    FacetGrid(std::shared_ptr<const DomainGridGeometry> gridGeometry,
              std::pair<std::unique_ptr<Grid>, std::vector<std::size_t>>&& pair)
    : domain_{std::move(gridGeometry)}
    , facetGrid_{std::move(pair.first)}
    , facetToDomainVertexIndex_{std::move(pair.second)}
    {}

    std::shared_ptr<const DomainGridGeometry> domain_;
    std::unique_ptr<Grid> facetGrid_;
    std::vector<std::size_t> facetToDomainVertexIndex_;
};

//! Convenience factory function for facet grids
template<typename Grid, typename GridGeometry, FacetSelectorFor<GridGeometry> Selector>
auto makeFacetGrid(std::shared_ptr<GridGeometry> gridGeometry, Selector&& selector)
{ return FacetGrid<Grid, std::remove_const_t<GridGeometry>>{std::move(gridGeometry), std::forward<Selector>(selector)}; }

/*!
 * \ingroup Discretization
 * \brief Holds a facet grid + connectivity mapping for finite-volume discretizations.
 * \tparam Grid The facet grid type
 * \tparam GridGeometry The grid geometry on which the facet grid is defined.
 */
template<typename Grid, typename GridGeometry>
class FVFacetGrid
{
    using FacetGrid = Dumux::FacetGrid<Grid, GridGeometry>;
    using FacetElementToScvfElementIndices = std::unordered_map<std::size_t, std::vector<std::size_t>>;
    static constexpr int domainDim = GridGeometry::GridView::dimension;

 public:
    using GridView = typename Grid::LeafGridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;
    using EntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    using DomainGridGeometry = GridGeometry;
    using DomainElement = typename GridGeometry::GridView::template Codim<0>::Entity;
    using DomainSubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    template<FacetSelectorFor<GridGeometry> FacetSelector>
    explicit FVFacetGrid(std::shared_ptr<const GridGeometry> gridGeometry, FacetSelector&& selector)
    : FVFacetGrid(std::make_shared<FacetGrid>(makeFacetGrid<Grid>(std::move(gridGeometry), std::forward<FacetSelector>(selector))))
    {}

    //! Constructor from a shared facet grid
    explicit FVFacetGrid(std::shared_ptr<const FacetGrid> facetGrid)
    : facetGrid_{std::move(facetGrid)}
    {
        const auto bboxTree = FacetGridDetail::makeBoundingBoxTree(gridView());

        domainElementToCouplingData_.resize(domainGridGeometry()->gridView().size(0));
        for (const auto& element : elements(domainGridGeometry()->gridView()))
        {
            // TODO: filter non-candidate elements to speed up computations?
            const auto eIdx = domainGridGeometry()->elementMapper().index(element);
            const auto fvGeometry = localView(*domainGridGeometry()).bindElement(element);
            for (const auto& scvf : scvfs(fvGeometry))
                if (
                    auto faceIdx = FacetGridDetail::coupledElementIndex(fvGeometry.geometry(scvf), gridView(), bboxTree);
                    faceIdx.has_value()
                )
                    domainElementToCouplingData_[eIdx][faceIdx.value()].push_back(scvf.index());
        }

        facetElementMapper_ = std::make_unique<EntityMapper>(gridView(), Dune::mcmgElementLayout());
        facetVertexMapper_ = std::make_unique<EntityMapper>(gridView(), Dune::mcmgVertexLayout());

        facetToDomainElements_.resize(gridView().size(0));
        for (std::size_t eIdxDomain = 0; eIdxDomain < domainGridGeometry()->gridView().size(0); ++eIdxDomain)
            for (const auto& [eIdxFacet, _] : domainElementToCouplingData_.at(eIdxDomain))
            {
                if (facetToDomainElements_[eIdxFacet].size() == 2)
                    DUNE_THROW(Dune::InvalidStateException, "Found more than two neighbors to a facet element");
                facetToDomainElements_[eIdxFacet].push_back(eIdxDomain);
            }
    }

    //! Return the grid view containing the selected facets
    GridView gridView() const
    { return facetGrid_->gridView(); }

    //! Return the index mapper for facet grid elements
    const EntityMapper& elementMapper() const
    { return *facetElementMapper_; }

    //! Return the index mapper for facet grid vertices
    const EntityMapper& vertexMapper() const
    { return *facetVertexMapper_; }

    //! Return the domain from which this facet grid was extracted
    std::shared_ptr<const GridGeometry> domainGridGeometry() const
    { return facetGrid_->domainGridGeometry(); }

    //! Return the index of the given vertex within the domain
    std::size_t domainVertexIndexOf(const Vertex& v) const
    { return facetGrid_->domainVertexIndexOf(v); }

    //! Return a range over all domain elements that overlap with the given facet grid element
    std::ranges::view auto domainElementsAdjacentTo(const Element& element) const
    {
        return facetToDomainElements_.at(elementMapper().index(element))
            | std::views::transform([&] (const auto& eIdxDomain) {
                return domainGridGeometry()->element(eIdxDomain);
            });
    }

    //! Return a range over the indices of the scvfs that overlap with the given trace element from within the given domain element
    std::ranges::view auto domainScvfsAdjacentTo(const Element& element, const DomainElement& domainElement) const
    {
        const auto eIdx = elementMapper().index(element);
        return domainElementToCouplingData_.at(domainGridGeometry()->elementMapper().index(domainElement))
            | std::views::filter([e=eIdx] (const auto& facetElementToScvfs) { return facetElementToScvfs.first == e; })
            | std::views::transform([&] (const auto& facetElementToScvfs) { return facetElementToScvfs.second; })
            | std::views::join;
    }

 private:
    std::shared_ptr<const FacetGrid> facetGrid_;
    std::unique_ptr<EntityMapper> facetElementMapper_;
    std::unique_ptr<EntityMapper> facetVertexMapper_;

    std::vector<FacetElementToScvfElementIndices> domainElementToCouplingData_;
    std::vector<Dune::ReservedVector<std::size_t, 2>> facetToDomainElements_;
};

//! Convenience factory function for finite-volume facet grids
template<typename Grid, typename GridGeometry, FacetSelectorFor<GridGeometry> Selector>
auto makeFVFacetGrid(std::shared_ptr<GridGeometry> gridGeometry, Selector&& selector)
{ return FVFacetGrid<Grid, std::remove_const_t<GridGeometry>>{std::move(gridGeometry), std::forward<Selector>(selector)}; }

/*!
 * \ingroup Discretization
 * \brief Extract the trace of a finite-volume discretization and exposes it as a new grid (or a subset of the trace).
 * \tparam GridGeometry The grid geometry from which to extract the facet grid.
 */
template<typename Grid, typename GridGeometry>
class FVTraceGrid : public FVFacetGrid<Grid, GridGeometry>
{
    using ParentType = FVFacetGrid<Grid, GridGeometry>;

 public:
    using typename ParentType::GridView;
    using typename ParentType::EntityMapper;

    using typename ParentType::DomainElement;
    using typename ParentType::DomainSubControlVolumeFace;

    //! Constructor for selecting a subset of the boundary
    template<FacetSelectorFor<GridGeometry> FacetSelector>
    explicit FVTraceGrid(std::shared_ptr<const GridGeometry> gg, FacetSelector&& selector)
    : ParentType{std::move(gg), [&] (const auto& is) { return is.boundary() and selector(is); }}
    {}

    //! Constructor for the trace of the entire boundary
    explicit FVTraceGrid(std::shared_ptr<const GridGeometry> gg)
    : ParentType{std::move(gg), [&] (const auto& is) { return is.boundary(); }}
    {}
};

//! Convenience factory function for finite-volume trace grids
template<typename Grid, typename GridGeometry, FacetSelectorFor<GridGeometry> Selector>
auto makeFVTraceGrid(std::shared_ptr<GridGeometry> gridGeometry, Selector&& selector)
{ return FVTraceGrid<Grid, std::remove_const_t<GridGeometry>>{std::move(gridGeometry), std::forward<Selector>(selector)}; }

//! Convenience factory function for finite-volume trace grids selecting the entire boundary
template<typename Grid, typename GridGeometry>
auto makeFVTraceGrid(std::shared_ptr<GridGeometry> gridGeometry)
{ return FVTraceGrid<Grid, std::remove_const_t<GridGeometry>>{std::move(gridGeometry)}; }

} // end namespace Dumux

#endif
