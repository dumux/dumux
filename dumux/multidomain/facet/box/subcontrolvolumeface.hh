// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
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
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/box/boxgeometryhelper.hh>
#include <dumux/discretization/box/subcontrolvolumeface.hh>

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
: public SubControlVolumeFaceBase<BoxFacetCouplingSubControlVolumeFace<GV, T>, T>
{
    using ThisType = BoxSubControlVolumeFace<GV, T>;
    using ParentType = SubControlVolumeFaceBase<ThisType, T>;
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
                                         const typename Element::Geometry& elemGeometry,
                                         unsigned int scvfIndex,
                                         std::vector<LocalIndexType>&& scvIndices)
    : corners_(geometryHelper.getScvfCorners(scvfIndex))
    , center_(0.0)
    , unitOuterNormal_(geometryHelper.normal(corners_, scvIndices))
    , area_(geometryHelper.scvfArea(corners_))
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , facetIndex_(/*undefined*/)
    , indexInFacet_(/*undefined*/)
    , boundary_(false)
    , interiorBoundary_(false)
    , boundaryFlag_{}
    {
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
    }

    //! Constructor for domain or interior boundary scvfs
    template<class GeometryHelper, class Intersection>
    BoxFacetCouplingSubControlVolumeFace(const GeometryHelper& geometryHelper,
                                         const Intersection& intersection,
                                         const typename Intersection::Geometry& isGeometry,
                                         LocalIndexType indexInIntersection,
                                         GridIndexType scvfIndex,
                                         std::vector<LocalIndexType>&& scvIndices,
                                         bool boundary,
                                         bool interiorBoundary)
    : corners_(geometryHelper.getBoundaryScvfCorners(intersection, isGeometry, indexInIntersection))
    , center_(0.0)
    , unitOuterNormal_(intersection.centerUnitOuterNormal())
    , area_(geometryHelper.scvfArea(corners_))
    , scvfIndex_(scvfIndex)
    , scvIndices_(std::move(scvIndices))
    , facetIndex_(intersection.indexInInside())
    , indexInFacet_(indexInIntersection)
    , boundary_(boundary)
    , interiorBoundary_(interiorBoundary)
    , boundaryFlag_{intersection}
    {
        for (const auto& corner : corners_)
            center_ += corner;
        center_ /= corners_.size();
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

    //! returns bolean if the sub control volume face is on the boundary
    bool boundary() const
    { return boundary_; }

    //! returns bolean if the sub control volume face is on an interior boundary
    bool interiorBoundary() const
    { return interiorBoundary_; }

    //! returns the unit nurmal vector pointing outwards
    //! of the sub-control volume that this scvf encloses
    const GlobalPosition& unitOuterNormal() const
    { return unitOuterNormal_; }

    //! index of the inside sub control volume for spatial param evaluation
    LocalIndexType insideScvIdx() const
    { return scvIndices_[0]; }

    //! The element-local index of this sub control volume face
    GridIndexType index() const
    { return scvfIndex_; }

    //! index of the outside sub control volume for spatial param evaluation
    //! This results in undefined behaviour if boundary is true
    LocalIndexType outsideScvIdx() const
    {
        assert(!boundary());
        return scvIndices_[1];
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

    //! The geometry of the sub control volume face
    Geometry geometry() const
    { return Geometry(Dune::GeometryTypes::cube(Geometry::mydimension), corners_); }

    //! Return the boundary flag
    typename BoundaryFlag::value_type boundaryFlag() const
    { return boundaryFlag_.get(); }

    //! returns the position of a corner of the face
    const GlobalPosition& corner(unsigned int localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

private:
    // geometrical information
    CornerStorage corners_;
    GlobalPosition center_;
    GlobalPosition unitOuterNormal_;
    Scalar area_;

    // indices
    GridIndexType scvfIndex_;
    std::vector<LocalIndexType> scvIndices_;

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
