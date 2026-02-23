// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1BubbleDiscretization
 * \brief Base class for the local finite element geometry for the pq1bubble method
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_FE_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_FE_ELEMENT_GEOMETRY_HH

#include <optional>
#include <utility>
#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include <dumux/discretization/pq1bubble/geometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Base class for the finite element geometry for pq1bubble discretizations.
 * \tparam GG the finite element grid geometry type
 * \tparam enableGridGeometryCache if the grid geometry is cached or not
 */
//! We don't need to do any specialization for caching since the local dof ranges are built on the fly.
template<class GG, bool enableGridGeometryCache>
class PQ1BubbleFEElementGeometry
{
    using GridView = typename GG::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GGCache = typename GG::Cache;
    using GeometryHelper = typename GGCache::GeometryHelper;

public:
    //! export the element type
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of finite element grid geometry
    using GridGeometry = GG;
    //! the quadrature rule type for elements
    using ElementQuadratureRule = typename GG::ElementQuadratureRule;
    //! the quadrature rule type for intersections
    using IntersectionQuadratureRule = typename GG::IntersectionQuadratureRule;
    //! the maximum number of scvs per element (2^dim for cubes)
    static constexpr std::size_t maxNumElementDofs = GridGeometry::maxNumElementDofs;

    //! Constructor
    PQ1BubbleFEElementGeometry(const GGCache& ggCache)
    : ggCache_(&ggCache)
    {}

    //! iterate over dof indices that are treated as hybrid dofs using the finite element method
    friend inline auto nonCVLocalDofs(const PQ1BubbleFEElementGeometry& fvGeometry)
    {
        return Dune::transformedRangeView(
            Dune::range(fvGeometry.numLocalDofs()),
            [&](const auto i) { return CVFE::LocalDof
            {
                static_cast<LocalIndexType>(i),
                static_cast<GridIndexType>(GeometryHelper::dofIndex(fvGeometry.gridGeometry().dofMapper(), fvGeometry.element(),
                                                                    fvGeometry.feLocalCoefficients().localKey(i))),
                static_cast<GridIndexType>(fvGeometry.elementIndex())
            }; }
        );
    }

    //! an iterator over all local dofs
    friend inline auto localDofs(const PQ1BubbleFEElementGeometry& fvGeometry)
    {
        return nonCVLocalDofs(fvGeometry);
    }

    //! an iterator over all local dofs related to an intersection
    template<class Intersection>
    friend inline auto localDofs(const PQ1BubbleFEElementGeometry& fvGeometry, const Intersection& intersection)
    {
        return Dune::transformedRangeView(
            Dune::range(GeometryHelper::numLocalDofsIntersection(fvGeometry.element().type(), intersection.indexInInside())),
            [&](const auto i)
            {
                auto localDofIdx = GeometryHelper::localDofIndexIntersection(fvGeometry.element().type(), intersection.indexInInside(), i);
                return CVFE::LocalDof
                {
                    static_cast<LocalIndexType>(localDofIdx),
                    static_cast<GridIndexType>(GeometryHelper::dofIndex(fvGeometry.gridGeometry().dofMapper(), fvGeometry.element(),
                                                                        fvGeometry.feLocalCoefficients().localKey(localDofIdx))),
                    static_cast<GridIndexType>(fvGeometry.elementIndex())
                };
                }
        );
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return gridGeometry().feCache().get(element_->type()).localBasis();
    }

    //! Get the local finite element coefficients
    const auto& feLocalCoefficients() const
    {
        return gridGeometry().feCache().get(element_->type()).localCoefficients();
    }

    //! The total number of element-local dofs
    std::size_t numLocalDofs() const
    {
        return GeometryHelper::numElementDofs(element().type());
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    PQ1BubbleFEElementGeometry bind(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! this function is for compatibility reasons with cc methods
    //! The pq1bubble stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element) &
    { this->bindElement(element); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    PQ1BubbleFEElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! Binding of an element, has to be called before using the fegeometries
    void bindElement(const Element& element) &
    {
        element_ = element;
        // cache element index
        eIdx_ = gridGeometry().elementMapper().index(element);
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

    //! The grid geometry we are a restriction of
    const GridGeometry& gridGeometry() const
    { return ggCache_->gridGeometry(); }

    //! The bound element index
    std::size_t elementIndex() const
    { return eIdx_; }

    //! Interpolation point data for a localDof
    template<class LocalDof>
    friend inline auto ipData(const PQ1BubbleFEElementGeometry& feGeometry, const LocalDof& localDof)
    {
        const auto type = feGeometry.element().type();
        const auto& localKey = feGeometry.gridGeometry().feCache().get(type).localCoefficients().localKey(localDof.index());
        const auto& localPos = GeometryHelper::localDofPosition(type, localKey);

        return CVFE::LocalDofInterpolationPointData{ localPos, feGeometry.elementGeometry().global(localPos), localDof.index() };
    }

    //! Interpolation point data for a global position
    friend inline auto ipData(const PQ1BubbleFEElementGeometry& feGeometry, const typename Element::Geometry::GlobalCoordinate& globalPos)
    {
        // Create ipData that does not automatically calculate the local position but only if it is called
        return CVFE::InterpolationPointDataLocalMapping{
            [&] (const typename Element::Geometry::GlobalCoordinate& pos) { return feGeometry.elementGeometry().local(pos); },
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
