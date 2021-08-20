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
 * \ingroup FaceCenteredStaggeredDiscretization
 * \copydoc Dumux::FreeFlowStaggeredSubControlVolumeFace
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUMEFACE_HH

#include <array>
#include <utility>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>
#include <dumux/discretization/facecentered/staggered/normalaxis.hh>

namespace Dumux {

/*!
 * \ingroup FaceCenteredStaggeredDiscretization
 * \brief Default traits class to be used for the sub-control volume face
 *        for the face-centered staggered finite volume scheme
 * \tparam GridView the type of the grid view
 */
template<class GridView>
struct FaceCenteredDefaultScvfGeometryTraits
{
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dim = GridView::Grid::dimension;
    static constexpr int dimWorld = GridView::Grid::dimensionworld;
    using CornerStorage = std::array<GlobalPosition, (1<<(dim-1))>;
    using Geometry = Dune::AxisAlignedCubeGeometry<Scalar, dim-1, dimWorld>;
};

template<class GridView, class T = FaceCenteredDefaultScvfGeometryTraits<GridView>>
class FaceCenteredStaggeredSubControlVolumeFace
{
    using Geometry = typename T::Geometry;
    using GridIndexType = typename T::GridIndexType;
    using Scalar = typename T::Scalar;
    using Element = typename T::Element;
    using CornerStorage = typename T::CornerStorage;

    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;

    using ElementGeometry = typename Element::Geometry;
    using IntersectionGeometry = typename GridView::Intersection::Geometry;

public:
    //! state the traits public and thus export all types
    using Traits = T;

    using GlobalPosition = typename T::GlobalPosition;
    enum class FaceType {frontal, lateral};

    FaceCenteredStaggeredSubControlVolumeFace() = default;

    //! The constructor for frontal faces
    FaceCenteredStaggeredSubControlVolumeFace(const ElementGeometry& elementGeometry,
                                              const IntersectionGeometry& intersectionGeometry,
                                              const std::array<GridIndexType, 2> globalScvIndices,
                                              const SmallLocalIndexType localScvfIdx,
                                              const GridIndexType globalScvfIdx,
                                              const GlobalPosition& unitOuterNormal,
                                              const FaceType faceType,
                                              const bool boundary)
    : globalScvIndices_(globalScvIndices)
    , localScvfIdx_(localScvfIdx)
    , globalScvfIdx_(globalScvfIdx)
    , area_(intersectionGeometry.volume())
    , normalAxis_(Dumux::normalAxis(unitOuterNormal))
    , outerNormalSign_(sign(unitOuterNormal[normalAxis_]))
    , faceType_(faceType)
    , boundary_(boundary)
    {
        assert(faceType == FaceType::frontal);
        center_ = boundary ? intersectionGeometry.center() : elementGeometry.center();
        ipGlobal_ = center_;

        if (!boundary)
            outerNormalSign_ *= -1.0;

        // the corners (coincide with intersection corners for boundary scvfs)
        const auto frontalOffSet = boundary ? GlobalPosition(0.0)
                                    : intersectionGeometry.center() - elementGeometry.center();

        for (int i = 0; i < corners_.size(); ++i)
            corners_[i] = intersectionGeometry.corner(i) + frontalOffSet;
    }

    //! The constructor for lateral faces
    template<class LateralFacetGeometry>
    FaceCenteredStaggeredSubControlVolumeFace(const ElementGeometry& elementGeometry,
                                              const IntersectionGeometry& intersectionGeometry,
                                              const LateralFacetGeometry& lateralFacetGeometry,
                                              const std::array<GridIndexType, 2> globalScvIndices,
                                              const SmallLocalIndexType localScvfIdx,
                                              const GridIndexType globalScvfIdx,
                                              const GlobalPosition& unitOuterNormal,
                                              const FaceType faceType,
                                              const bool boundary)
    : globalScvIndices_(globalScvIndices)
    , localScvfIdx_(localScvfIdx)
    , globalScvfIdx_(globalScvfIdx)
    , area_(0.5*lateralFacetGeometry.volume())
    , normalAxis_(Dumux::normalAxis(unitOuterNormal))
    , outerNormalSign_(sign(unitOuterNormal[normalAxis_]))
    , faceType_(faceType)
    , boundary_(boundary)
    {
        assert(faceType == FaceType::lateral);
        const auto shift = intersectionGeometry.center() - elementGeometry.center();
        ipGlobal_ = lateralFacetGeometry.center() + shift;
        center_ = 0.5*(lateralFacetGeometry.center() + ipGlobal_);

        const auto dofAxis = Dumux::normalAxis(shift);
        const auto eps = shift.two_norm() * 1e-8;

        for (int i = 0; i < corners_.size(); ++i)
        {
            // copy the corner of the corresponding lateral facet
            auto& corner = corners_[i];
            corner = lateralFacetGeometry.corner(i);

            // shift the corner such that the scvf covers half of the lateral facet
            // (keep the outer corner positions)
            using std::abs;
            if (abs(corner[dofAxis] - intersectionGeometry.center()[dofAxis]) > eps)
                corner[dofAxis] = elementGeometry.center()[dofAxis];
        }
    }

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    { return center_; }

    //! The integration point of the sub control volume face
    const GlobalPosition& ipGlobal() const
    { return ipGlobal_; }

    //! The unit outer normal
    const GlobalPosition unitOuterNormal() const
    {
        GlobalPosition result(0.0);
        result[normalAxis_] = 1.0 * directionSign();
        return result;
    }

    //! Index of the inside sub control volume for spatial param evaluation
    GridIndexType insideScvIdx() const
    { return globalScvIndices_[0]; }

    //! index of the outside sub control volume for spatial param evaluation
    GridIndexType outsideScvIdx() const
    { return globalScvIndices_[1]; }

    GridIndexType index() const
    { return globalScvfIdx_; }

    SmallLocalIndexType localIndex() const
    { return localScvfIdx_; }

    FaceType faceType() const
    { return faceType_; }

    bool boundary() const
    { return boundary_; }

    bool isFrontal() const
    { return faceType_ == FaceType::frontal; }

    bool isLateral() const
    { return faceType_ == FaceType::lateral; }

    Scalar area() const
    { return area_; }

    SmallLocalIndexType normalAxis() const
    { return normalAxis_; }

    std::int_least8_t directionSign() const
    { return outerNormalSign_; }

    const GlobalPosition& corner(unsigned int localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

    //! The geometry of the sub control volume face
    Geometry geometry() const
    {
        auto inPlaneAxes = std::move(std::bitset<T::dimWorld>{}.set());
        inPlaneAxes.set(normalAxis_, false);
        return { corners_.front(), corners_.back(), inPlaneAxes };
    }

private:
    GlobalPosition center_;
    GlobalPosition ipGlobal_;
    CornerStorage corners_;
    std::array<GridIndexType, 2> globalScvIndices_;
    SmallLocalIndexType localScvfIdx_;
    GridIndexType globalScvfIdx_;
    Scalar area_;
    SmallLocalIndexType normalAxis_;
    std::int_least8_t outerNormalSign_;
    FaceType faceType_;
    bool boundary_;
};

} // end namespace Dumux

#endif
