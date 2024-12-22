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
#include <type_traits>
#include <unordered_map>
#include <algorithm>
#include <ranges>

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace FacetGridMapperDetail {

    template<typename Geometry, typename FacetGridView, typename FacetBoundingBoxTree>
    auto overlappingFacetElementIndices(const Geometry& geometry,
                                        const FacetGridView& facetGridView,
                                        const FacetBoundingBoxTree& bboxTree)
    {
        std::vector<std::size_t> result;
        const auto intersections = intersectingEntities(geometry, bboxTree);
        if (intersections.empty())
            return result;
        result.resize(intersections.size());
        std::ranges::copy(intersections | std::views::transform([&] (const auto& is) { return is.second(); }), result.begin());
        std::ranges::sort(result);
        result.erase(std::unique(result.begin(), result.end()), result.end());
        return result;
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

        for (const auto& element : elements(domainGridGeometry_->gridView()))
        {
            // TODO: filter non-candidate elements to speed up computations?
            const auto eIdx = domainGridGeometry_->elementMapper().index(element);
            const auto fvGeometry = localView(*domainGridGeometry_).bindElement(element);
            if constexpr (isCVFE)
            {
                for (const auto& scv : scvs(fvGeometry))
                    for (const auto facetElementIndex : FacetGridMapperDetail::overlappingFacetElementIndices(
                        fvGeometry.geometry(scv),
                        facetGridView,
                        bboxTree
                    ))
                        // TODO: localDofIndex robust?
                        domainElementToCouplingData_[eIdx][facetElementIndex].push_back(scv.localDofIndex());
            }
            else
            {
                for (const auto& scvf : scvfs(fvGeometry))
                    for (const auto facetElementIndex : FacetGridMapperDetail::overlappingFacetElementIndices(
                        fvGeometry.geometry(scvf),
                        facetGridView,
                        bboxTree
                    ))
                        domainElementToCouplingData_[eIdx][facetElementIndex].push_back(scvf.index());
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
        return domainScesAdjacentTo_(element, domainElement);
    }

    //! Return a range over the indices of the scvs that overlap with the given trace element from within the given domain element
    std::ranges::view auto domainScvsAdjacentTo(const FacetElement& element, const DomainElement& domainElement) const
    {
        static_assert(isCVFE, "Scv mapping only available for cvfe methods");
        return domainScesAdjacentTo_(element, domainElement);
    }

 private:
    std::ranges::view auto domainScesAdjacentTo_(const FacetElement& element, const DomainElement& domainElement) const
    {
        const auto eIdx = facetEntitySet_->index(element);
        return domainElementToCouplingData_.at(domainGridGeometry_->elementMapper().index(domainElement))
            | std::views::filter([e=eIdx] (const auto& facetElementToScvfs) { return facetElementToScvfs.first == e; })
            | std::views::transform([&] (const auto& facetElementToScvfs) { return facetElementToScvfs.second; })
            | std::views::join;
    }

    FacetGridView facetGridView_;
    std::shared_ptr<FacetEntitySet> facetEntitySet_;
    std::shared_ptr<const DomainGridGeometry> domainGridGeometry_;
    std::vector<FacetElementToScvfElementIndices> domainElementToCouplingData_;
    std::vector<Dune::ReservedVector<std::size_t, 2>> facetToDomainElements_;
};

template<typename FGV, typename GG>
FVFacetGridMapper(const FGV&, std::shared_ptr<GG>) -> FVFacetGridMapper<FGV, std::remove_const_t<GG>>;

} // end namespace Dumux

#endif
