// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PQ1BubbleDiscretization
 * \brief Base class for the local finite volume geometry for the pq1bubble method
 *        This builds up the sub control volumes and sub control volume faces
 *        for an element.
 */
#ifndef DUMUX_DISCRETIZATION_PQ1BUBBLE_FV_ELEMENT_GEOMETRY_HH
#define DUMUX_DISCRETIZATION_PQ1BUBBLE_FV_ELEMENT_GEOMETRY_HH

#include <optional>
#include <utility>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/pqkfactory.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/cvfe/integrationpointdata.hh>

#include <dumux/discretization/pq1bubble/geometryhelper.hh>

namespace Dumux {

/*!
 * \ingroup PQ1BubbleDiscretization
 * \brief Base class for the finite volume geometry vector for pq1bubble models
 *        This builds up the sub control volumes and sub control volume faces
 *        for each element.
 * \tparam GG the finite volume grid geometry type
 * \tparam enableGridGeometryCache if the grid geometry is cached or not
 */
template<class GG, bool enableGridGeometryCache>
class PQ1BubbleFVElementGeometry;

//! specialization in case the FVElementGeometries are stored
template<class GG>
class PQ1BubbleFVElementGeometry<GG, true>
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
    using IpData = Dumux::CVFE::LocalDofIntegrationPointData<typename GridView::template Codim<0>::Entity::Geometry::LocalCoordinate,
                                                             typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate,
                                                             LocalIndexType>;

public:
    //! export the element type
    using Element = typename GridView::template Codim<0>::Entity;
    //! export type of subcontrol volume
    using SubControlVolume = typename GG::SubControlVolume;
    //! export type of subcontrol volume face
    using SubControlVolumeFace = typename GG::SubControlVolumeFace;
    //! export type of finite volume grid geometry
    using GridGeometry = GG;
    //! the maximum number of scvs per element (2^dim for cubes)
    // ToDo get this from GG
    static constexpr std::size_t maxNumElementScvs = (1<<dim) + 1;

    //! Constructor
    PQ1BubbleFVElementGeometry(const GGCache& ggCache)
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
    scvs(const PQ1BubbleFVElementGeometry& fvGeometry)
    {
        using Iter = typename std::vector<SubControlVolume>::const_iterator;
        const auto& s = fvGeometry.ggCache_->scvs(fvGeometry.eIdx_);
        return Dune::IteratorRange<Iter>(s.begin(), s.end());
    }

    //! iterate over dof indices that belong to dofs associated with control volumes
    friend inline auto cvLocalDofs(const PQ1BubbleFVElementGeometry& fvGeometry)
    {
        return Dune::transformedRangeView(
            Dune::range(fvGeometry.numLocalDofs()-GeometryHelper::numNonCVLocalDofs(fvGeometry.element().type())),
            [&](const auto i) { return CVFE::LocalDof
            {
                static_cast<LocalIndexType>(i),
                static_cast<GridIndexType>(GeometryHelper::dofIndex(fvGeometry.gridGeometry().dofMapper(), fvGeometry.element(), i)),
                static_cast<GridIndexType>(fvGeometry.elementIndex())
            }; }
        );
    }

    //! iterate over dof indices that are treated as hybrid dofs using the finite element method
    template<bool enable = GridGeometry::enableHybridCVFE, std::enable_if_t<enable, int> = 0>
    friend inline auto nonCVLocalDofs(const PQ1BubbleFVElementGeometry& fvGeometry)
    {
        return Dune::transformedRangeView(
            Dune::range(fvGeometry.numLocalDofs()-GeometryHelper::numNonCVLocalDofs(fvGeometry.element().type()), fvGeometry.numLocalDofs()),
            [&](const auto i) { return CVFE::LocalDof
            {
                static_cast<LocalIndexType>(i),
                static_cast<GridIndexType>(GeometryHelper::dofIndex(fvGeometry.gridGeometry().dofMapper(), fvGeometry.element(), i)),
                static_cast<GridIndexType>(fvGeometry.elementIndex())
            }; }
        );
    }

    //! an iterator over all local dofs
    friend inline auto localDofs(const PQ1BubbleFVElementGeometry& fvGeometry)
    {
        return Dune::transformedRangeView(
            Dune::range(fvGeometry.numLocalDofs()),
            [&](const auto i) { return CVFE::LocalDof
            {
                static_cast<LocalIndexType>(i),
                static_cast<GridIndexType>(GeometryHelper::dofIndex(fvGeometry.gridGeometry().dofMapper(), fvGeometry.element(), i)),
                static_cast<GridIndexType>(fvGeometry.elementIndex())
            }; }
        );
    }

    //! an iterator over all local dofs related to an intersection
    template<class Intersection>
    friend inline auto localDofs(const PQ1BubbleFVElementGeometry& fvGeometry, const Intersection& intersection)
    {
        return Dune::transformedRangeView(
            Dune::range(GeometryHelper::numLocalDofsIntersection(fvGeometry.element().type(), intersection.indexInInside())),
            [&](const auto i)
            {
                auto localDofIdx = GeometryHelper::localDofIndexIntersection(fvGeometry.element().type(), intersection.indexInInside(), i);
                return CVFE::LocalDof
                {
                    static_cast<LocalIndexType>(localDofIdx),
                    static_cast<GridIndexType>(GeometryHelper::dofIndex(fvGeometry.gridGeometry().dofMapper(), fvGeometry.element(), localDofIdx)),
                    static_cast<GridIndexType>(fvGeometry.elementIndex())
                };
                }
        );
    }

    //! get the position related to a localdof
    template<class LocalDof>
    auto dofPosition(const LocalDof& localDof) const
    {
        return GeometryHelper::dofPosition(this->element(), localDof.index());
    }

    //! iterator range for sub control volumes faces. Iterates over
    //! all scvfs of the bound element.
    //! This is a free function found by means of ADL
    //! To iterate over all sub control volume faces of this FVElementGeometry use
    //! for (auto&& scvf : scvfs(fvGeometry))
    friend inline Dune::IteratorRange<typename std::vector<SubControlVolumeFace>::const_iterator>
    scvfs(const PQ1BubbleFVElementGeometry& fvGeometry)
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

    //! The total number of element-local dofs
    std::size_t numLocalDofs() const
    {
        return GeometryHelper::numElementDofs(element().type());
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
    PQ1BubbleFVElementGeometry bind(const Element& element) &&
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
    PQ1BubbleFVElementGeometry bindElement(const Element& element) &&
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
    }

    //! Returns true if bind/bindElement has already been called
    bool isBound() const
    { return static_cast<bool>(element_); }

    //! The bound element
    const Element& element() const
    { return *element_; }

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
        if (scv.isOverlapping())
            DUNE_THROW(Dune::NotImplemented, "Geometry of overlapping scv");

        assert(isBound());
        const auto geo = element().geometry();
        const GeometryHelper helper(geo);
        return {
            helper.getScvGeometryType(scv.indexInElement()),
            helper.getScvCorners(scv.indexInElement())
        };
    }

    //! Geometry of a sub control volume face
    typename SubControlVolumeFace::Traits::Geometry geometry(const SubControlVolumeFace& scvf) const
    {
        assert(isBound());
        const auto geo = element().geometry();
        if (scvf.boundary())
        {
            GeometryHelper helper(geo);
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
            GeometryHelper helper(geo);
            return {
                helper.getInteriorScvfGeometryType(scvf.index()),
                helper.getScvfCorners(scvf.index())
            };
        }
    }

    //! Integration point data for an scv
    friend inline IpData ipData(const PQ1BubbleFVElementGeometry& fvGeometry, const SubControlVolume& scv)
    {
        const auto type = fvGeometry.element().type();
        const auto& localKey = fvGeometry.gridGeometry().feCache().get(type).localCoefficients().localKey(scv.localDofIndex());

        return IpData(GeometryHelper::localDofPosition(type, localKey), scv.dofPosition(), scv.localDofIndex());
    }

private:
    const GGCache* ggCache_;
    GridIndexType eIdx_;

    std::optional<Element> element_;
};

} // end namespace Dumux

#endif
