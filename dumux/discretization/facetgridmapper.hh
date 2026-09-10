// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \brief Mapping between the elements of a facet grid and the intersections of the domain grid
 */
#ifndef DUMUX_DISCRETIZATION_FACET_GRID_MAPPER_HH
#define DUMUX_DISCRETIZATION_FACET_GRID_MAPPER_HH

#include <limits>
#include <memory>
#include <ranges>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/reservedvector.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/localview.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Maps between the elements of a grid defined on the facets of a domain grid and the
 *        intersections of the domain grid they were extracted from.
 *
 * A facet grid element is an intersection of the domain grid, with a domain element on each
 * side. A discretization exposes such an intersection on the domain boundary as a boundary
 * face of its local view, and everything defined on it, the sub-control volume faces, the
 * local degrees of freedom and the quadrature rule, is obtained from the local view.
 *
 * \tparam FacetGridView The facet grid view type
 * \tparam GG The grid geometry of the domain
 */
template<typename FacetGridView, typename GG>
class FacetGridMapper
{
    using GridIndexType = typename IndexTraits<typename GG::GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<typename GG::GridView>::LocalIndex;

    static constexpr bool isCVFE = DiscretizationMethods::isCVFE<typename GG::DiscretizationMethod>;
    static constexpr int domainDim = GG::GridView::dimension;
    static constexpr int facetDim = domainDim - 1;
    static_assert(int(FacetGridView::dimension) == facetDim);
    static_assert(int(FacetGridView::dimensionworld) == GG::GridView::dimensionworld);

 public:
    using DomainGridGeometry = GG;
    using DomainElement = typename DomainGridGeometry::GridView::template Codim<0>::Entity;
    using FacetElement = typename FacetGridView::template Codim<0>::Entity;
    using FacetVertex = typename FacetGridView::template Codim<facetDim>::Entity;
    using BoundaryFace = typename DomainGridGeometry::BoundaryFace;

    //! A side of a facet element: a domain element and the index of the intersection within it
    struct Side
    {
        GridIndexType elementIndex;
        LocalIndexType intersectionIndex;
    };

    template<typename FacetGridManager>
    FacetGridMapper(const FacetGridManager& facetGridManager, std::shared_ptr<const DomainGridGeometry> gridGeometry)
    : facetGridView_{facetGridManager.grid().leafGridView()}
    , domainGridGeometry_{std::move(gridGeometry)}
    {
        sides_.resize(facetGridView_.size(0));
        for (const auto& facetElement : elements(facetGridView_))
            for (const auto& record : facetGridManager.hostGridIntersections(facetElement))
                sides_[facetGridView_.indexSet().index(facetElement)].push_back({
                    static_cast<GridIndexType>(record.elementIndex),
                    static_cast<LocalIndexType>(record.indexInInside)
                });
    }

    //! Return a range over all domain elements adjacent to the given facet grid element
    std::ranges::view auto domainElementsAdjacentTo(const FacetElement& element) const
    {
        return sides_[facetGridView_.indexSet().index(element)]
            | std::views::transform([&] (const Side& side) { return domainGridGeometry_->element(side.elementIndex); });
    }

    //! Return the index, within the given domain element, of the intersection the facet element lies on
    LocalIndexType intersectionIndex(const FacetElement& element, const DomainElement& domainElement) const
    {
        const auto eIdx = domainGridGeometry_->elementMapper().index(domainElement);
        for (const auto& side : sides_[facetGridView_.indexSet().index(element)])
            if (side.elementIndex == eIdx)
                return side.intersectionIndex;
        DUNE_THROW(Dune::InvalidStateException, "The given domain element is not adjacent to the facet element");
    }

    /*!
     * \brief Return the boundary face of the given local view that the facet element lies on.
     * \note The local view must be bound to a domain element adjacent to the facet element, and
     *       the facet element must lie on the boundary of the domain.
     */
    template<typename FVElementGeometry>
    BoundaryFace boundaryFace(const FVElementGeometry& fvGeometry, const FacetElement& element) const
    {
        const auto isIdx = intersectionIndex(element, fvGeometry.element());
        for (const auto& face : boundaryFaces(fvGeometry))
            if (face.intersectionIndex() == isIdx)
                return face;
        DUNE_THROW(Dune::InvalidStateException, "The facet element does not lie on the boundary of the domain");
    }

    //! Return the indices of the sub-control volume faces on the given facet element within the given domain element
    std::vector<GridIndexType> domainScvfsAdjacentTo(const FacetElement& element, const DomainElement& domainElement) const
    {
        const auto fvGeometry = localView(*domainGridGeometry_).bindElement(domainElement);
        std::vector<GridIndexType> result;
        // the faces of a control-volume finite element scheme lie inside the elements, so only
        // a boundary face carries any; a cell-centered scheme has faces on every intersection
        if constexpr (isCVFE)
        {
            const auto isIdx = intersectionIndex(element, domainElement);
            for (const auto& face : boundaryFaces(fvGeometry))
                if (face.intersectionIndex() == isIdx)
                    for (const auto& scvf : scvfs(fvGeometry, face))
                        result.push_back(scvf.index());
        }
        else
            for (const auto& scvf : scvfs(fvGeometry, faceOf_(fvGeometry, element)))
                result.push_back(scvf.index());
        return result;
    }

    /*!
     * \brief Return the indices of the local degrees of freedom on the given facet element within
     *        the given domain element, or the degree of freedom of the element for schemes whose
     *        degrees of freedom are not associated with the boundary.
     */
    std::vector<std::size_t> domainLocalDofsAdjacentTo(const FacetElement& element, const DomainElement& domainElement) const
    {
        const auto fvGeometry = localView(*domainGridGeometry_).bindElement(domainElement);
        std::vector<std::size_t> result;
        if constexpr (requires { localDofs(fvGeometry, std::declval<const BoundaryFace&>()); })
            for (const auto& localDof : localDofs(fvGeometry, faceOf_(fvGeometry, element)))
                result.push_back(localDof.index());
        else
            for (const auto& scv : scvs(fvGeometry))
                result.push_back(scv.dofIndex());
        return result;
    }

    //! Return the grid geometry of the domain
    const DomainGridGeometry& domainGridGeometry() const
    { return *domainGridGeometry_; }

 private:
    // The boundary face the facet element lies on or, for an interior facet, a face built from
    // the intersection itself. The latter carries no local index and serves only queries that
    // are answered from the intersection index.
    template<typename FVElementGeometry>
    BoundaryFace faceOf_(const FVElementGeometry& fvGeometry, const FacetElement& element) const
    {
        const auto isIdx = intersectionIndex(element, fvGeometry.element());
        for (const auto& face : boundaryFaces(fvGeometry))
            if (face.intersectionIndex() == isIdx)
                return face;

        for (const auto& is : intersections(domainGridGeometry_->gridView(), fvGeometry.element()))
            if (is.indexInInside() == isIdx)
                return BoundaryFace{
                    is.geometry().center(),
                    is.geometry().volume(),
                    is.centerUnitOuterNormal(),
                    std::numeric_limits<LocalIndexType>::max(),
                    isIdx,
                    typename BoundaryFace::Traits::BoundaryFlag{is}
                };
        DUNE_THROW(Dune::InvalidStateException, "The domain element has no intersection with the given index");
    }

    FacetGridView facetGridView_;
    std::shared_ptr<const DomainGridGeometry> domainGridGeometry_;
    std::vector<Dune::ReservedVector<Side, 2>> sides_;
};

template<typename FacetGridManager, typename GG>
FacetGridMapper(const FacetGridManager&, std::shared_ptr<GG>)
    -> FacetGridMapper<typename FacetGridManager::Grid::LeafGridView, std::remove_const_t<GG>>;

} // end namespace Dumux

#endif
