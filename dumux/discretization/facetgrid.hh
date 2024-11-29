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
#include <ranges>

#include <dune/common/reservedvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/mcmgmapper.hh>

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif  // HAVE_DUNE_FOAMGRID

#include <dumux/common/typetraits/typetraits.hh>
#include <dumux/common/typetraits/problem.hh>

#include <dumux/discretization/extrusion.hh>
#include <dumux/geometry/geometryintersection.hh>

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
    struct _DefaultFacetGrid
    {
#if HAVE_DUNE_FOAMGRID
        using type = Dune::FoamGrid<GridView::dimension - 1, GridView::dimension>;
#else
        static_assert(AlwaysFalse<GridView>::value, "FoamGrid not found; please explicitly specify a GridType for the trace.");
#endif
    };

    template<typename GridGeometry>
    using DefaultFacetGrid = typename _DefaultFacetGrid<typename GridGeometry::GridView>::type;

}  // namespace FacetGridDetail
#endif  // DOXYGEN

template<typename T, typename GridGeometry>
concept FacetSelectorFor = std::invocable<T, const typename GridGeometry::GridView::Intersection&>
    and std::is_convertible_v<std::invoke_result_t<T, const typename GridGeometry::GridView::Intersection&>, bool>;

/*!
 * \ingroup Discretization
 * \brief Extracts facets of a finite-volume discretization and exposes them as a new grid.
 * \tparam GridGeometry The grid geometry from which to extract the facet grid.
 */
template<typename GridGeometry, typename FacetGrid = FacetGridDetail::DefaultFacetGrid<GridGeometry>>
class FVFacetGrid
{
    using FacetElementToScvfElementIndices = std::unordered_map<std::size_t, std::vector<std::size_t>>;
    static constexpr int domainDim = GridGeometry::GridView::dimension;

 public:
    using GridView = typename FacetGrid::LeafGridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using EntityMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;

    using DomainElement = typename GridGeometry::GridView::template Codim<0>::Entity;
    using DomainSubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;

    template<FacetSelectorFor<GridGeometry> FacetSelector>
    explicit FVFacetGrid(std::shared_ptr<const GridGeometry> gg, FacetSelector&& selector)
    : domainGridGeometry_{std::move(gg)}
    { update(std::forward<FacetSelector>(selector)); }

    //! Recompute the facet grid and mappings (e.g. after grid adaptation or with different selector)
    template<FacetSelectorFor<GridGeometry> FacetSelector>
    void update(FacetSelector&& selector)
    {
        facetGridFactory_ = std::make_unique<Dune::GridFactory<FacetGrid>>();

        std::vector<std::vector<unsigned int>> elementCorners;
        std::vector<unsigned int> localCornerStorage;

        elementCorners.reserve(domainGridView_().size(0));
        domainElementToCouplingData_.resize(domainGridView_().size(0));
        auto fvGeometry = localView(*domainGridGeometry_);

        std::ranges::for_each(elements(domainGridView_()), [&] (const auto& element) {
            const auto eIdx = domainGridGeometry_->elementMapper().index(element);
            const auto& refElement = Dune::referenceElement(element);
            const auto& elemGeo = element.geometry();

            auto getFVGeometry = [&, isBound=false] () mutable -> const typename GridGeometry::LocalView& {
                if (!isBound) { fvGeometry.bindElement(element); isBound = true; }
                return fvGeometry;
            };

            std::ranges::for_each(
                intersections(domainGridView_(), element) | std::views::filter(selector),
                [&] (const auto& intersection) {
                    const auto fvGeometry = getFVGeometry();
                    const auto& isGeo = intersection.geometry();

                    bool mayBeDuplicateElement = true;
                    localCornerStorage.clear();
                    localCornerStorage.reserve(isGeo.corners());
                    std::ranges::for_each(std::views::iota(0, isGeo.corners()), [&] (int c) {
                        const auto vIdxLocal = refElement.subEntity(intersection.indexInInside(), 1, c, domainDim);
                        const auto vIdxGlobal = domainGridGeometry_->vertexMapper().subIndex(element, vIdxLocal, domainDim);
                        localCornerStorage.push_back([&] () {
                            auto it = domainToFacetVertexMap_.find(vIdxGlobal);
                            if (it == domainToFacetVertexMap_.end())
                            {
                                mayBeDuplicateElement = false;
                                const auto currentVertexIndex = domainToFacetVertexMap_.size();
                                facetGridFactory_->insertVertex(elemGeo.global(refElement.position(vIdxLocal, domainDim)));
                                it = domainToFacetVertexMap_.emplace(static_cast<std::size_t>(vIdxGlobal), currentVertexIndex).first;
                            }
                            return it->second;
                        } ());
                    });

                    std::optional<std::size_t> facetElementIndex;
                    if (mayBeDuplicateElement) {
                        auto it = std::find_if(elementCorners.begin(), elementCorners.end(), [&] (const auto& corners) {
                            const auto contained = [&] (auto idx) -> bool { return std::ranges::count(corners, idx); };
                            return std::ranges::all_of(localCornerStorage, contained);
                        });
                        if (it != elementCorners.end())
                            facetElementIndex = std::distance(elementCorners.begin(), it);
                    }

                    if (!facetElementIndex.has_value())
                    {
                        facetElementIndex = elementCorners.size();
                        elementCorners.push_back(localCornerStorage);
                        facetGridFactory_->insertElement(isGeo.type(), localCornerStorage);
                    }

                    std::ranges::for_each(scvfs(fvGeometry), [&] (const auto& scvf) {
                        if (FacetGridDetail::intersect(fvGeometry.geometry(scvf), isGeo))
                            domainElementToCouplingData_[eIdx][facetElementIndex.value()].push_back(scvf.index());
                    });
            });
        });

        facetGrid_ = facetGridFactory_->createGrid();
        facetElementMapper_ = std::make_unique<EntityMapper>(gridView(), Dune::mcmgElementLayout());
        facetVertexMapper_ = std::make_unique<EntityMapper>(gridView(), Dune::mcmgVertexLayout());

        // map insertion to mapper indices for scvf map
        std::vector<std::size_t> insertionToActualIndex(facetGrid_->leafGridView().size(0));
        for (const auto& e : elements(facetGrid_->leafGridView()))
            insertionToActualIndex[facetGridFactory_->insertionIndex(e)] = facetElementMapper_->index(e);
        for (auto& map : domainElementToCouplingData_) {
            FacetElementToScvfElementIndices new_map;
            for (const auto& [facetElemInsertionIndex, scvfIndices] : map)
                new_map[insertionToActualIndex[facetElemInsertionIndex]] = scvfIndices;
            map = std::move(new_map);
        }

        facetToDomainElements_.resize(gridView().size(0));
        for (std::size_t eIdxDomain = 0; eIdxDomain < domainGridView_().size(0); ++eIdxDomain)
            for (const auto& [eIdxFacet, _] : domainElementToCouplingData_.at(eIdxDomain))
                facetToDomainElements_[eIdxFacet].push_back(eIdxDomain);
    }

    //! Return the grid view containing the selected facets
    GridView gridView() const
    { return facetGrid_->leafGridView(); }

    //! Return the index mapper for facet grid elements
    const EntityMapper& elementMapper() const
    { return *facetElementMapper_; }

    //! Return the index mapper for facet grid vertices
    const EntityMapper& vertexMapper() const
    { return *facetVertexMapper_; }

    //! Return a range over all domain elements that overlap with the given facet grid element
    std::ranges::view auto adjacentDomainElements(const Element& element) const
    {
        return facetToDomainElements_.at(elementMapper().index(element))
            | std::views::transform([&] (const auto& eIdxDomain) {
                return domainGridGeometry_->element(eIdxDomain);
            });
    }

    //! Return a range over the indices of the scvfs that overlap with the given trace element from within the given domain element
    std::ranges::view auto adjacentScvfIndices(const Element& element, const DomainElement& domainElement) const
    {
        const auto eIdx = elementMapper().index(element);
        return domainElementToCouplingData_.at(domainGridGeometry_->elementMapper().index(domainElement))
            | std::views::filter([e=eIdx] (const auto& facetElementToScvfs) { return facetElementToScvfs.first == e; })
            | std::views::transform([&] (const auto& facetElementToScvfs) { return facetElementToScvfs.second; })
            | std::views::join;
    }

    //! Return the index of vertex that coincides with the given domain vertex
    std::size_t domainToFacetVertexIndex(std::size_t domainVertexIndex) const
    { return domainToFacetVertexMap_.at(domainVertexIndex); }

    //! Return a reference to the grid geometry that this grad has been extracted from
    const GridGeometry& domainGridGeometry() const
    { return *domainGridGeometry_; }

 private:
    const auto& domainGridView_() const { return domainGridGeometry_->gridView(); }

    std::shared_ptr<const GridGeometry> domainGridGeometry_;

    std::unique_ptr<Dune::GridFactory<FacetGrid>> facetGridFactory_;
    std::unique_ptr<FacetGrid> facetGrid_;

    std::unique_ptr<EntityMapper> facetElementMapper_;
    std::unique_ptr<EntityMapper> facetVertexMapper_;

    std::vector<FacetElementToScvfElementIndices> domainElementToCouplingData_;
    std::vector<Dune::ReservedVector<std::size_t, 2>> facetToDomainElements_;
    std::unordered_map<std::size_t, std::size_t> domainToFacetVertexMap_;
};

template<typename GridGeometry, typename Indicator>
FVFacetGrid(std::shared_ptr<GridGeometry>, Indicator&&) -> FVFacetGrid<std::remove_const_t<GridGeometry>>;


/*!
 * \ingroup Discretization
 * \brief Extract the trace of a finite-volume discretization and exposes it as a new grid (or a subset of the trace).
 * \tparam GridGeometry The grid geometry from which to extract the facet grid.
 */
template<typename GridGeometry, typename FacetGrid = FacetGridDetail::DefaultFacetGrid<GridGeometry>>
class FVTraceGrid : private FVFacetGrid<GridGeometry, FacetGrid>
{
    using ParentType = FVFacetGrid<GridGeometry, FacetGrid>;

 public:
    using typename ParentType::GridView;
    using typename ParentType::EntityMapper;

    using typename ParentType::DomainElement;
    using typename ParentType::DomainSubControlVolumeFace;

    template<FacetSelectorFor<GridGeometry> FacetSelector>
    explicit FVTraceGrid(std::shared_ptr<const GridGeometry> gg, FacetSelector&& selector)
    : ParentType{std::move(gg), [&] (const auto& is) { return is.boundary() and selector(is); }}
    {}

    //! Constructor
    explicit FVTraceGrid(std::shared_ptr<const GridGeometry> gg)
    : ParentType{std::move(gg), [&] (const auto& is) { return is.boundary(); }}
    {}

    //! Recompute the facet grid and mappings (e.g. after grid adaptation or with different selector)
    template<FacetSelectorFor<GridGeometry> FacetSelector>
    void update(FacetSelector&& selector)
    { ParentType::update([&] (const auto& is) { return is.boundary() and selector(is); }); }

    //! Recompute the facet grid and mappings (e.g. after grid adaptation)
    void update()
    { ParentType::update([&] (const auto& is) { return is.boundary(); }); }

    using ParentType::gridView;
    using ParentType::elementMapper;
    using ParentType::vertexMapper;
    using ParentType::domainGridGeometry;
    using ParentType::adjacentDomainElements;
    using ParentType::adjacentScvfIndices;
    using ParentType::domainToFacetVertexIndex;
};

template<typename GridGeometry>
FVTraceGrid(std::shared_ptr<GridGeometry>) -> FVTraceGrid<std::remove_const_t<GridGeometry>>;

template<typename GridGeometry, typename Indicator>
FVTraceGrid(std::shared_ptr<GridGeometry>, Indicator&&) -> FVTraceGrid<std::remove_const_t<GridGeometry>>;

} // end namespace Dumux

#endif
