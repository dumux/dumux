// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BoxDiscretization
 * \brief Base class for the local finite element discretization for the PQ1 FE method.
 *        Treats all vertex degrees of freedom as FE dofs.
 */
#ifndef DUMUX_DISCRETIZATION_PQ1_FE_ELEMENT_DISCRETIZATION_HH
#define DUMUX_DISCRETIZATION_PQ1_FE_ELEMENT_DISCRETIZATION_HH

#include <optional>
#include <ranges>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/common/rangeutilities.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include <dumux/discretization/pq2/geometryhelper.hh> // for FEPQ2GeometryHelper (reused for P1)

namespace Dumux {

/*!
 * \ingroup BoxDiscretization
 * \brief Local finite element discretization for the PQ1 FE method.
 *        All degrees of freedom (vertices only) are treated as FE dofs.
 * \tparam GG the finite element grid discretization type (PQ1FEGridDiscretization)
 * \tparam enableGridDiscretizationCache if the grid discretization is cached or not
 */
template<class GG, bool enableGridDiscretizationCache>
class PQ1FEElementDiscretization
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GGCache = typename GG::Cache;
    using GeometryHelper = typename GGCache::GeometryHelper;

public:
    //! export the element type
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of finite element grid discretization
    using GridDiscretization = GG;
    using GridGeometry [[deprecated("Use GridDiscretization instead")]] = GG;
    //! the quadrature rule type for elements
    using ElementQuadratureRule = typename GG::ElementQuadratureRule;
    //! the quadrature rule type for intersections
    using IntersectionQuadratureRule = typename GG::IntersectionQuadratureRule;
    //! the quadrature rule type for boundary faces
    using BoundaryFaceQuadratureRule = typename GG::BoundaryFaceQuadratureRule;
    //! export the boundary face type
    using BoundaryFace = typename GG::BoundaryFace;
    //! the maximum number of dofs per element (2^dim for cubes)
    static constexpr std::size_t maxNumElementDofs = GridDiscretization::maxNumElementDofs;

    //! Constructor
    PQ1FEElementDiscretization(const GGCache& ggCache)
    : ggCache_(&ggCache)
    {}

    //! iterate over all local dofs (vertex shape function dofs)
    friend inline auto localDofs(const PQ1FEElementDiscretization& elemDisc)
    {
        return Dune::transformedRangeView(
            Dune::range(elemDisc.numLocalDofs()),
            [&](const auto i) {
                return CVFE::LocalDof{
                    static_cast<LocalIndexType>(i),
                    static_cast<GridIndexType>(GeometryHelper::dofIndex(
                        elemDisc.gridDiscretization().dofMapper(),
                        elemDisc.element(),
                        elemDisc.feLocalCoefficients().localKey(i))),
                    static_cast<GridIndexType>(elemDisc.elementIndex())
                };
            }
        );
    }

    //! iterate over non-CV local dofs (for the pure FE model all dofs are non-CV)
    friend inline auto nonCVLocalDofs(const PQ1FEElementDiscretization& elemDisc)
    {
        return localDofs(elemDisc);
    }

    //! an iterator over all local dofs related to a boundary face
    friend inline auto localDofs(const PQ1FEElementDiscretization& elemDisc, const BoundaryFace& boundaryFace)
    {
        return std::views::iota(std::size_t(0), elemDisc.numLocalDofs())
            | std::views::filter([&](std::size_t i) {
                return GeometryHelper::localDofOnIntersection(
                    elemDisc.element().type(),
                    boundaryFace.intersectionIndex(),
                    elemDisc.feLocalCoefficients().localKey(i));
            })
            | std::views::transform([&](std::size_t i) {
                return CVFE::LocalDof{
                    static_cast<LocalIndexType>(i),
                    static_cast<GridIndexType>(GeometryHelper::dofIndex(
                        elemDisc.gridDiscretization().dofMapper(),
                        elemDisc.element(),
                        elemDisc.feLocalCoefficients().localKey(i))),
                    static_cast<GridIndexType>(elemDisc.elementIndex())
                };
            });
    }

    //! iterator range for boundary faces of the bound element.
    //! To iterate: for (auto&& bf : boundaryFaces(elemDisc))
    friend inline std::ranges::view auto
    boundaryFaces(const PQ1FEElementDiscretization& elemDisc)
    {
        const auto& v = elemDisc.ggCache_->boundaryFaces(elemDisc.eIdx_);
        return std::ranges::views::all(v);
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return gridDiscretization().feCache().get(element_->type()).localBasis();
    }

    //! Get the local finite element coefficients
    const auto& feLocalCoefficients() const
    {
        return gridDiscretization().feCache().get(element_->type()).localCoefficients();
    }

    //! The total number of element-local dofs
    std::size_t numLocalDofs() const
    {
        return feLocalCoefficients().size();
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    PQ1FEElementDiscretization bind(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! this function is for compatibility reasons with cc methods
    //! The pq1 stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element) &
    { this->bindElement(element); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    PQ1FEElementDiscretization bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! Binding of an element, has to be called before using the fe geometries
    void bindElement(const Element& element) &
    {
        element_ = element;
        eIdx_ = gridDiscretization().elementMapper().index(element);
        elementGeometry_.emplace(element.geometry());
    }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

    //! The bound element geometry
    const typename Element::Geometry& elementGeometry() const
    { return *elementGeometry_; }

    //! The grid geometry we are a restriction of (deprecated)
    [[deprecated("Use gridDiscretization() instead")]]
    const GridDiscretization& gridGeometry() const
    { return ggCache_->gridDiscretization(); }

    //! The grid discretization we are a restriction of
    const GridDiscretization& gridDiscretization() const
    { return ggCache_->gridDiscretization(); }

    //! Returns whether the element has boundary faces
    bool hasBoundaryFaces() const
    { return !ggCache_->boundaryFaces(eIdx_).empty(); }

    //! Get a boundary face with a local boundary face index
    const BoundaryFace& boundaryFace(LocalIndexType bfIdx) const
    { return ggCache_->boundaryFaces(eIdx_)[bfIdx]; }

    //! The bound element index
    std::size_t elementIndex() const
    { return eIdx_; }

    //! Geometry of a boundary face
    typename BoundaryFace::Traits::Geometry geometry(const BoundaryFace& boundaryFace) const
    {
        assert(isBound());
        const auto& elemGeo = elementGeometry();
        const auto faceGeoInRef = referenceElement(elemGeo).template geometry<1>(boundaryFace.intersectionIndex());
        typename BoundaryFace::Traits::CornerStorage corners;
        for (int i = 0; i < faceGeoInRef.corners(); ++i)
            corners.push_back(elemGeo.global(faceGeoInRef.corner(i)));
        return { faceGeoInRef.type(), corners };
    }

    //! Interpolation point data for a localDof
    template<class LocalDof>
    friend inline auto ipData(const PQ1FEElementDiscretization& elemDisc, const LocalDof& localDof)
    {
        const auto type = elemDisc.element().type();
        const auto& localKey = elemDisc.feLocalCoefficients().localKey(localDof.index());
        const auto& localPos = GeometryHelper::localDofPosition(type, localKey);
        return CVFE::LocalDofInterpolationPointData{
            localPos,
            elemDisc.elementGeometry().global(localPos),
            localDof.index()
        };
    }

    //! Interpolation point data for a global position
    friend inline auto ipData(const PQ1FEElementDiscretization& elemDisc,
                              const typename Element::Geometry::GlobalCoordinate& globalPos)
    {
        return CVFE::InterpolationPointDataLocalMapping{
            [&] (const typename Element::Geometry::GlobalCoordinate& pos) {
                return elemDisc.elementGeometry().local(pos); },
            globalPos
        };
    }

private:
    const GGCache* ggCache_;
    GridIndexType eIdx_;

    std::optional<Element> element_;
    std::optional<typename Element::Geometry> elementGeometry_;
};

} // end namespace Dumux

#endif
