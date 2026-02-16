// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ2Discretization
 * \brief Base class for the local finite volume geometry for the pq2 method
 *        This builds up the sub control volumes and sub control volume faces
 *        for an element.
 */
#ifndef DUMUX_DISCRETIZATION_PQ2_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_PQ2_FV_ELEMENT_GEOMETRY_HH

#include <cstddef>
#include <cassert>
#include <optional>
#include <ranges>
#include <utility>
#include <vector>

#include <dune/common/rangeutilities.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/interpolationpointdata.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>

#include <dumux/discretization/pq2/geometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup PQ2Discretization
 * \brief Base class for the finite volume geometry vector for pq2 models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 * \tparam GG the finite volume grid geometry type
 * \tparam enableGridGeometryCache if the grid geometry is cached or not
 */
template<class GG, bool enableGridGeometryCache>
class PQ2FVElementGeometry;

//! specialization in case the FVElementGeometries are stored
template<class GG>
class PQ2FVElementGeometry<GG, true>
{
    using GridView = typename GG::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using CoordScalar = typename GridView::ctype;
    using FeLocalBasis = typename GG::FeCache::FiniteElementType::Traits::LocalBasisType;
    using GGCache = typename GG::Cache;
    using GeometryHelper = typename GGCache::GeometryHelper;

    using BaseIpData = CVFE::InterpolationPointData<
                        typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate,
                        typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate
                        >;

public:
    //! export the element type
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    //! the quadrature rule type for scvs
    using ScvQuadratureRule = typename GG::ScvQuadratureRule;
    //! the quadrature rule type for scvfs
    using ScvfQuadratureRule = typename GG::ScvfQuadratureRule;
    //! the quadrature rule type for elements
    using ElementQuadratureRule = typename GG::ElementQuadratureRule;
    //! the quadrature rule type for intersections
    using IntersectionQuadratureRule = typename GG::IntersectionQuadratureRule;
    //! the maximum number of scvs per element for hypercubes
    static constexpr std::size_t maxNumElementDofs = GridGeometry::maxNumElementDofs;

    //! Constructor
    PQ2FVElementGeometry(const GGCache& ggCache)
    : ggCache_(&ggCache)
    {}

    //! Get a sub control volume with a local scv index
    const SubControlVolume& scv(LocalIndexType scvIdx) const
    {
        return ggCache_->scvs(eIdx_)[scvIdx];
    }

    //! Get a sub control volume face with a local scvf index
    const SubControlVolumeFace& scvf(LocalIndexType scvfIdx) const
    {
        return ggCache_->scvfs(eIdx_)[scvfIdx];
    }

    //! iterator range for sub control volumes. Iterates over
    //! all scvs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volumes of this FVElementGeometry use
    //! for (auto&& scv : scvs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolume>::const_iterator>
    scvs(const PQ2FVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        const auto& s = fvGeometry.ggCache_->scvs(fvGeometry.eIdx_);
        return Dune::IteratorRange<Iter>(s.begin(), s.end());
    }

    friend inline auto cvLocalDofs(const PQ2FVElementGeometry& fvGeometry)
    {
        return std::views::iota(std::size_t(0), fvGeometry.numLocalDofs())
            | std::views::filter([&](size_t i) { return fvGeometry.feLocalCoefficients().localKey(i).codim() == dim; })
            | std::views::transform([&](size_t i) {
                return CVFE::LocalDof{
                    static_cast<LocalIndexType>(i),
                    static_cast<GridIndexType>(GeometryHelper::dofIndex(fvGeometry.gridGeometry().dofMapper(), fvGeometry.element(),
                                                                        fvGeometry.feLocalCoefficients().localKey(i))),
                    static_cast<GridIndexType>(fvGeometry.elementIndex())
                };
            });
    }

    friend inline auto nonCVLocalDofs(const PQ2FVElementGeometry& fvGeometry)
    {
        return std::views::iota(std::size_t(0), fvGeometry.numLocalDofs())
            | std::views::filter([&](size_t i) { return !(fvGeometry.feLocalCoefficients().localKey(i).codim() == dim); })
            | std::views::transform([&](size_t i) {
                return CVFE::LocalDof{
                    static_cast<LocalIndexType>(i),
                    static_cast<GridIndexType>(GeometryHelper::dofIndex(fvGeometry.gridGeometry().dofMapper(), fvGeometry.element(),
                                                                        fvGeometry.feLocalCoefficients().localKey(i))),
                    static_cast<GridIndexType>(fvGeometry.elementIndex())
                };
            });
    }

    friend inline auto localDofs(const PQ2FVElementGeometry& fvGeometry)
    {
        return Dune::transformedRangeView(
            Dune::range(std::size_t(0), fvGeometry.numLocalDofs()), [&](const auto i) {
                return CVFE::LocalDof{
                    static_cast<LocalIndexType>(i),
                    static_cast<GridIndexType>(GeometryHelper::dofIndex(fvGeometry.gridGeometry().dofMapper(), fvGeometry.element(),
                                                                        fvGeometry.feLocalCoefficients().localKey(i))),
                    static_cast<GridIndexType>(fvGeometry.elementIndex())
                };
            }
        );
    }

    //! an iterator over all local dofs related to an intersection
    template<class Intersection>
    friend inline auto localDofs(const PQ2FVElementGeometry& fvGeometry, const Intersection& intersection)
    {
        return std::views::iota(std::size_t(0), fvGeometry.numLocalDofs())
            | std::views::filter([&](size_t i) {
                return GeometryHelper::localDofOnIntersection(fvGeometry.element().type(),
                                                              intersection.indexInInside(),
                                                              fvGeometry.feLocalCoefficients().localKey(i)); })
            | std::views::transform([&](size_t i) {
                return CVFE::LocalDof{
                    static_cast<LocalIndexType>(i),
                    static_cast<GridIndexType>(GeometryHelper::dofIndex(fvGeometry.gridGeometry().dofMapper(), fvGeometry.element(),
                                                                        fvGeometry.feLocalCoefficients().localKey(i))),
                    static_cast<GridIndexType>(fvGeometry.elementIndex())
                };
            });
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const PQ2FVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolumeFace>::const_iterator;
        const auto& s = fvGeometry.ggCache_->scvfs(fvGeometry.eIdx_);
        return Dune::IteratorRange<Iter>(s.begin(), s.end());
    }

    //! Get a local finite element basis
    const FeLocalBasis& feLocalBasis() const
    {
        return gridGeometry().feCache().get(element_->type()).localBasis();
    }

    //! Get a local finite element basis
    const auto& feLocalCoefficients() const
    {
        return gridGeometry().feCache().get(element_->type()).localCoefficients();
    }

    //! The total number of element-local dofs
    std::size_t numLocalDofs() const
    {
        return feLocalCoefficients().size();
    }

    //! The total number of sub control volumes
    std::size_t numScv() const
    {
        return ggCache_->scvs(eIdx_).size();
    }

    //! The total number of sub control volume faces
    std::size_t numScvf() const
    {
        return ggCache_->scvfs(eIdx_).size();
    }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bind(element);`
     */
    PQ2FVElementGeometry bind(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! this function is for compatibility reasons with cc methods
    //! The pq2 stencil is always element-local so bind and bindElement
    //! are identical.
    void bind(const Element& element) &
    { this->bindElement(element); }

    /*!
     * \brief bind the local view (r-value overload)
     * This overload is called when an instance of this class is a temporary in the usage context
     * This allows a usage like this: `const auto view = localView(...).bindElement(element);`
     */
    PQ2FVElementGeometry bindElement(const Element& element) &&
    {
        this->bindElement(element);
        return std::move(*this);
    }

    //! Binding of an element, has to be called before using the fvgeometries
    //! Prepares all the volume variables within the element
    //! For compatibility reasons with the FVGeometry cache being disabled
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

    //! Returns whether one of the geometry's scvfs lies on a boundary
    bool hasBoundaryScvf() const
    { return ggCache_->hasBoundaryScvf(eIdx_); }

    //! The bound element index
    std::size_t elementIndex() const
    { return eIdx_; }

    //! The intersection index the scvf belongs to
    std::size_t intersectionIndex(const SubControlVolumeFace& scvf) const
    {
        const auto localScvfIdx = scvf.index() - GeometryHelper::numInteriorScvf(element().type());
        return ggCache_->scvfBoundaryGeometryKeys(eIdx_)[localScvfIdx][0];
    }

    //! Geometry of a sub control volume
    typename SubControlVolume::Traits::Geometry geometry(const SubControlVolume& scv) const
    {
        assert(isBound());
        const auto& geo = elementGeometry();
        const GeometryHelper helper(geo);
        const auto scvIdx = this->feLocalCoefficients().localKey(scv.localDofIndex()).subEntity();
        return {
            helper.getScvGeometryType(scvIdx),
            helper.getScvCorners(scvIdx)
        };
    }

    //! Geometry of a sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    {
        assert(isBound());
        if (scvf.boundary())
        {
            GeometryHelper helper(*elementGeometry_);
            const auto localScvfIdx = scvf.index() - GeometryHelper::numInteriorScvf(element().type());
            const auto [localFacetIndex, isScvfLocalIdx]
                = ggCache_->scvfBoundaryGeometryKeys(eIdx_)[localScvfIdx];
            return {
                helper.getBoundaryScvfGeometryType(isScvfLocalIdx),
                helper.getBoundaryScvfCorners(localFacetIndex, isScvfLocalIdx)
            };
        }
        else
        {
            GeometryHelper helper(*elementGeometry_);
            return {
                helper.getInteriorScvfGeometryType(scvf.index()),
                helper.getScvfCorners(scvf.index())
            };
        }
    }

    //! Interpolation point data for an scv
    friend inline auto ipData(const PQ2FVElementGeometry& fvGeometry, const SubControlVolume& scv)
    {
        const auto type = fvGeometry.element().type();
        const auto& localKey = fvGeometry.feLocalCoefficients().localKey(scv.localDofIndex());

        return CVFE::LocalDofInterpolationPointData{ GeometryHelper::localDofPosition(type, localKey), scv.dofPosition(), scv.localDofIndex() };
    }

    //! Interpolation point data for a localDof
    template<class LocalDof>
    friend inline auto ipData(const PQ2FVElementGeometry& fvGeometry, const LocalDof& localDof)
    {
        const auto type = fvGeometry.element().type();
        const auto& localKey = fvGeometry.feLocalCoefficients().localKey(localDof.index());
        const auto& localPos = GeometryHelper::localDofPosition(type, localKey);

        return CVFE::LocalDofInterpolationPointData{ localPos, fvGeometry.elementGeometry().global(localPos), localDof.index() };
    }

    //! Interpolation point data for a global position
    friend inline auto ipData(const PQ2FVElementGeometry& fvGeometry, const typename Element::Geometry::GlobalCoordinate& globalPos)
    {
        // Create ipData that does not automatically calculate the local position but only if it is called
        return CVFE::InterpolationPointDataLocalMapping{
            [&] (const typename Element::Geometry::GlobalCoordinate& pos) { return fvGeometry.elementGeometry().local(pos); },
            globalPos
        };
    }

    //! Interpolation point data for scvf
    friend inline auto ipData(const PQ2FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf)
    {
        return CVFE::FaceInterpolationPointData<BaseIpData, LocalIndexType>
                    { scvf.unitOuterNormal(), scvf.index(), fvGeometry.elementGeometry().local(scvf.ipGlobal()), scvf.ipGlobal() };
    }

private:
    const GGCache* ggCache_;
    GridIndexType eIdx_;

    std::optional<Element> element_;
    std::optional<typename Element::Geometry> elementGeometry_;
};

} // end namespace Dumux

#endif
