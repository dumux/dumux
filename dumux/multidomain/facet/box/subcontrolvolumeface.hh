// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetCoupling
 * \brief Base class for a sub control volume face of the box method
 *        in the context of of models considering coupling of different
 *        domains across the bulk grid facets.
 */
#ifndef DUMUX_FACETCOUPLING_BOX_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_FACETCOUPLING_BOX_SUBCONTROLVOLUMEFACE_HH

#include <utility>

#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/boundaryflag.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/discretization/box/subcontrolvolumeface.hh>
#include <dumux/geometry/volume.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Class for a sub control volume face in the box method, i.e a part of the boundary
 *        of a sub control volume we compute fluxes on. This is a specialization for models
 *        considering coupling of different domains across the bulk grid facets.
 * \tparam GV the type of the grid view
 * \tparam T the scvf geometry traits
 */
template<class GV, class T = BoxDefaultScvfGeometryTraits<GV> >
class BoxFacetCouplingSubControlVolumeFace
{
    using GridIndexType = typename T::GridIndexType;
    using LocalIndexType = typename T::LocalIndexType;
    using Scalar = typename T::Scalar;
    using CornerStorage = typename T::CornerStorage;
    using Geometry = typename T::Geometry;
    using BoundaryFlag = typename T::BoundaryFlag;

public:
    //! export the type used for global coordinates
    using GlobalPosition = typename T::GlobalPosition;
    //! state the traits public and thus export all types
    using Traits = T;

    //! The default constructor
    BoxFacetCouplingSubControlVolumeFace() = default;

    //! Constructor for inner scvfs
    template<class GeometryHelper, class Element>
    BoxFacetCouplingSubControlVolumeFace(const GeometryHelper& geometryHelper,
                                         const Element& element,
                                         unsigned int scvfIndex,
                                         std::array<LocalIndexType, 2>&& scvIndices)
    : center_(0.0)
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , facetIndex_(/*undefined*/)
    , indexInFacet_(/*undefined*/)
    , boundary_(false)
    , interiorBoundary_(false)
    , boundaryFlag_{}
    {
        const auto corners = geometryHelper.getScvfCorners(scvfIndex);
        unitOuterNormal_ = geometryHelper.normal(corners, scvIndices_);
        area_ = Dumux::convexPolytopeVolume<T::dim-1>(
            Dune::GeometryTypes::cube(T::dim-1),
            [&](unsigned int i){ return corners[i]; });
        for (const auto& corner : corners)
            center_ += corner;
        center_ /= corners.size();
    }

    //! Constructor for domain or interior boundary scvfs
    template<class GeometryHelper, class Intersection>
    BoxFacetCouplingSubControlVolumeFace(const GeometryHelper& geometryHelper,
                                         const Intersection& intersection,
                                         LocalIndexType indexInIntersection,
                                         GridIndexType scvfIndex,
                                         std::array<LocalIndexType, 2>&& scvIndices,
                                         bool boundary,
                                         bool interiorBoundary)
    : center_(0.0)
    , unitOuterNormal_(intersection.centerUnitOuterNormal())
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , facetIndex_(intersection.indexInInside())
    , indexInFacet_(indexInIntersection)
    , boundary_(boundary)
    , interiorBoundary_(interiorBoundary)
    , boundaryFlag_{intersection}
    {
        auto corners = geometryHelper.getBoundaryScvfCorners(intersection.indexInInside(), indexInIntersection);
        area_ = Dumux::convexPolytopeVolume<T::dim-1>(
            Dune::GeometryTypes::cube(T::dim-1),
            [&](unsigned int i){ return corners[i]; });
        for (const auto& corner : corners)
            center_ += corner;
        center_ /= corners.size();
    }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    { return center_; }

    //! The integration point for flux evaluations in global coordinates
    const GlobalPosition& ipGlobal() const
    { return center_; }

    //! The area of the sub control volume face
    Scalar area() const
    { return area_; }

    //! returns true if the sub control volume face is on the boundary
    bool boundary() const
    { return boundary_; }

    //! returns true if the sub control volume face is on an interior boundary
    bool interiorBoundary() const
    { return interiorBoundary_; }

    //! returns the unit nurmal vector pointing outwards
    //! of the sub-control volume that this scvf encloses
    const GlobalPosition& unitOuterNormal() const
    { return unitOuterNormal_; }

    //! index of the inside sub control volume
    LocalIndexType insideScvIdx() const
    { return scvIndices_[0]; }

    //! The element-local index of this sub control volume face
    GridIndexType index() const
    { return scvfIndex_; }

    //! Index of the i-th outside sub control volume or boundary scv index.
    // Results in undefined behaviour if i >= numOutsideScvs()
    LocalIndexType outsideScvIdx(int i = 0) const
    {
        assert(!boundary() && !interiorBoundary());
        return scvIndices_[1];
    }

    //! The number of scvs on the outside of this face
    std::size_t numOutsideScvs() const
    {
        return static_cast<std::size_t>(!(boundary() || interiorBoundary()));
    }

    //! returns the element-local index of the facet this scvf is embedded in.
    //! This is only valid to be called for scvfs on domain/interior boundaries.
    LocalIndexType facetIndexInElement() const
    {
        assert(interiorBoundary_ || boundary_);
        return facetIndex_;
    }

    //! returns the facet-local (intersection-local) index of this scvf.
    //! This is only valid to be called for scvfs on domain/interior boundaries.
    LocalIndexType indexInElementFacet() const
    {
        assert(interiorBoundary_ || boundary_);
        return indexInFacet_;
    }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    { return boundaryFlag_.get(); }

private:
    // geometrical information
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;

    // indices
    GridIndexType scvfIndex_;
    std::array<LocalIndexType, 2> scvIndices_;

    // indices valid for domain/interior boundary scvfs
    LocalIndexType facetIndex_;
    LocalIndexType indexInFacet_;

    // boundary information
    bool boundary_;
    bool interiorBoundary_;
    BoundaryFlag boundaryFlag_;
};

} // end namespace Dumux

#endif
