// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Maps between entities of a discretization and grids defined on its facets.
 */
#ifndef DUMUX_DISCRETIZATION_FACET_GRID_MAPPER_HH
#define DUMUX_DISCRETIZATION_FACET_GRID_MAPPER_HH

#include <vector>
#include <memory>
#include <utility>
#include <optional>
#include <type_traits>
#include <unordered_map>
#include <ranges>

#include <dune/common/reservedvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/geometry/geometryintersection.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace FacetGridMapperDetail {

    template<typename Geo1, typename Geo2>
    bool intersect(const Geo1& geo1, const Geo2& geo2)
    {
        using Algo = GeometryIntersection<Geo1, Geo2>;
        typename Algo::Intersection r;
        return Algo::intersection(geo1, geo2, r);
    }

    template<typename ScvfGeo, typename FacetGridView, typename FacetBoundingBoxTree>
    auto overlappingFacetElement(const ScvfGeo& scvfGeo,
                                 const FacetGridView& facetGridView,
                                 const FacetBoundingBoxTree& bboxTree)
    {
        using FacetElement = typename FacetGridView::template Codim<0>::Entity;
        // TODO: use `intersectingEntities` once it has been made robust for bboxes with width=0 in one direction
        // TODO: assert that only one -> only conforming grids implemented at this point
        for (const auto& element : elements(facetGridView))
            if (intersect(scvfGeo, element.geometry()))
                return std::optional<FacetElement>{std::in_place, element};
        return std::optional<FacetElement>{};
    }

}  // namespace FacetGridMapperDetail
#endif  // DOXYGEN

/*!
 * \ingroup Discretization
 * \brief Maps between entities of finite-volume discretizations and
 *        a grid defined on the facets of the discretization.
 * \tparam FacetGridView The facet grid view type
 * \tparam GridGeometry The grid geometry on which the facet grid is defined.
 */
template<typename FacetGridView, typename GG>
class FVFacetGridMapper
{
    using FacetEntitySet = GridViewGeometricEntitySet<FacetGridView>;
    using FacetElementToScvfElementIndices = std::unordered_map<std::size_t, std::vector<std::size_t>>;
    static constexpr int domainDim = GG::GridView::dimension;
    static constexpr int facetDim = domainDim - 1;
    static constexpr bool isCVFE = Dumux::DiscretizationMethods::isCVFE<typename GG::DiscretizationMethod>;
    static_assert(int(FacetGridView::dimension) == facetDim);
    static_assert(int(FacetGridView::dimensionworld) == GG::GridView::dimensionworld);

 public:
    using DomainGridGeometry = GG;
    using DomainElement = typename DomainGridGeometry::GridView::template Codim<0>::Entity;
    using FacetElement = typename FacetGridView::template Codim<0>::Entity;
    using FacetVertex = typename FacetGridView::template Codim<facetDim>::Entity;

    explicit FVFacetGridMapper(const FacetGridView& facetGridView, std::shared_ptr<const DomainGridGeometry> gridGeometry)
    : facetGridView_{facetGridView}
    , facetEntitySet_{std::make_shared<FacetEntitySet>(facetGridView)}
    , domainGridGeometry_{std::move(gridGeometry)}
    {
        BoundingBoxTree<FacetEntitySet> bboxTree{facetEntitySet_};
        domainElementToCouplingData_.resize(domainGridGeometry_->gridView().size(0));

        if (isCVFE)
            facetToDomainVertex_.resize(facetGridView.size(FacetGridView::dimension), std::size_t{0});

        for (const auto& element : elements(domainGridGeometry_->gridView()))
        {
            // TODO: filter non-candidate elements to speed up computations?
            const auto eIdx = domainGridGeometry_->elementMapper().index(element);
            const auto fvGeometry = localView(*domainGridGeometry_).bindElement(element);
            if constexpr (isCVFE)
            {
                for (const auto& scv : scvs(fvGeometry))
                    if (
                        auto facetElement = FacetGridMapperDetail::overlappingFacetElement(
                            fvGeometry.geometry(scv),
                            facetGridView,
                            bboxTree
                        );
                        facetElement.has_value()
                    )
                    {
                        // TODO: localDofIndex robust?
                        domainElementToCouplingData_[eIdx][facetEntitySet_->index(*facetElement)].push_back(scv.localDofIndex());
                        const auto& facetElementGeo = facetElement->geometry();
                        const auto& facetRefElement = Dune::referenceElement(*facetElement);
                        const auto& dofPosition = scv.dofPosition();
                        for (unsigned int corner = 0; corner < facetElement->subEntities(facetDim); ++corner)
                            if (Dune::FloatCmp::eq(
                                facetElementGeo.global(facetRefElement.position(corner, facetDim)),
                                dofPosition
                            ))
                                // TODO: `dofIndex` robust??
                                facetToDomainVertex_[facetGridView.indexSet().subIndex(*facetElement, corner, facetDim)] = scv.dofIndex();
                    }
            }
            else
            {
                for (const auto& scvf : scvfs(fvGeometry))
                    if (
                        auto facetElement = FacetGridMapperDetail::overlappingFacetElement(
                            fvGeometry.geometry(scvf),
                            facetGridView,
                            bboxTree
                        );
                        facetElement.has_value()
                    )
                        domainElementToCouplingData_[eIdx][facetEntitySet_->index(*facetElement)].push_back(scvf.index());
            }
        }

        facetToDomainElements_.resize(domainGridGeometry_->gridView().size(0));
        for (std::size_t eIdxDomain = 0; eIdxDomain < domainGridGeometry_->gridView().size(0); ++eIdxDomain)
            for (const auto& [eIdxFacet, _] : domainElementToCouplingData_.at(eIdxDomain))
            {
                if (facetToDomainElements_[eIdxFacet].size() == 2)
                    DUNE_THROW(Dune::InvalidStateException, "Found more than two neighbors to a facet element");
                facetToDomainElements_[eIdxFacet].push_back(eIdxDomain);
            }
    }

    //! Return the index of the given vertex within the domain
    std::size_t domainVertexIndexOf(const FacetVertex& v) const
    {
        static_assert(isCVFE, "Vertex mapping currently only implemented for CVFE methods");
        return facetToDomainVertex_.at(facetGridView_.indexSet().index(v));
    }

    //! Return a range over all domain elements that overlap with the given facet grid element
    std::ranges::view auto domainElementsAdjacentTo(const FacetElement& element) const
    {
        return facetToDomainElements_.at(facetEntitySet_->index(element))
            | std::views::transform([&] (const auto& eIdxDomain) {
                return domainGridGeometry_->element(eIdxDomain);
            });
    }

    //! Return a range over the indices of the scvfs that overlap with the given trace element from within the given domain element
    std::ranges::view auto domainScvfsAdjacentTo(const FacetElement& element, const DomainElement& domainElement) const
    {
        static_assert(!isCVFE, "Scvf mapping unavailable for cvfe methods");
        const auto eIdx = facetEntitySet_->index(element);
        return domainElementToCouplingData_.at(domainGridGeometry_->elementMapper().index(domainElement))
            | std::views::filter([e=eIdx] (const auto& facetElementToScvfs) { return facetElementToScvfs.first == e; })
            | std::views::transform([&] (const auto& facetElementToScvfs) { return facetElementToScvfs.second; })
            | std::views::join;
    }

    //! Return a range over the indices of the scvs that overlap with the given trace element from within the given domain element
    std::ranges::view auto domainScvsAdjacentTo(const FacetElement& element, const DomainElement& domainElement) const
    {
        static_assert(isCVFE, "Scv mapping only available for cvfe methods");
        const auto eIdx = facetEntitySet_->index(element);
        return domainElementToCouplingData_.at(domainGridGeometry_->elementMapper().index(domainElement))
            | std::views::filter([e=eIdx] (const auto& facetElementToScvfs) { return facetElementToScvfs.first == e; })
            | std::views::transform([&] (const auto& facetElementToScvfs) { return facetElementToScvfs.second; })
            | std::views::join;
    }

 private:
    FacetGridView facetGridView_;
    std::shared_ptr<FacetEntitySet> facetEntitySet_;
    std::shared_ptr<const DomainGridGeometry> domainGridGeometry_;
    std::vector<FacetElementToScvfElementIndices> domainElementToCouplingData_;
    std::vector<Dune::ReservedVector<std::size_t, 2>> facetToDomainElements_;
    std::vector<std::size_t> facetToDomainVertex_;
};

template<typename FGV, typename GG>
FVFacetGridMapper(const FGV&, std::shared_ptr<GG>) -> FVFacetGridMapper<FGV, std::remove_const_t<GG>>;

} // end namespace Dumux

#endif
