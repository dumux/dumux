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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::FreeFlowStaggeredSubControlVolumeFace
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUMEFACE_HH

#include <array>
#include <utility>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/discretization/subcontrolvolumefacebase.hh>


#include <typeinfo>

namespace Dumux {

template <class GridView>
class FaceCenteredStaggeredSubControlVolumeFace
{
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GridView::ctype;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;

public:
    using GlobalPosition =typename Element::Geometry::GlobalCoordinate;
    enum class FaceType {frontal, lateral};

    FaceCenteredStaggeredSubControlVolumeFace() = default;

    FaceCenteredStaggeredSubControlVolumeFace(const GlobalPosition& center,
                                              const GlobalPosition& ipGlobal,
                                              const std::array<SmallLocalIndexType, 2> localScvIndices,
                                              const std::array<GridIndexType, 2> globalScvIndices,
                                              const SmallLocalIndexType localScvfIdx,
                                              const Scalar area,
                                              const SmallLocalIndexType directionIdx,
                                              const std::int_least8_t outerNormalSign,
                                              const GridIndexType globalScvfIdx,
                                              const GridIndexType scvfIdxWithCommonEntity,
                                              const FaceType faceType,
                                              const bool boundary)
    : center_(center)
    , ipGlobal_(ipGlobal)
    , localScvIndices_(localScvIndices)
    , globalScvIndices_(globalScvIndices)
    , localScvfIdx_(localScvfIdx)
    , area_(area)
    , directionIdx_(directionIdx)
    , outerNormalSign_(outerNormalSign)
    , globalScvfIdx_(globalScvfIdx)
    , scvfIdxWithCommonEntity_(scvfIdxWithCommonEntity)
    , faceType_(faceType)
    , boundary_(boundary) {}

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
        result[directionIndex()] = 1.0 * directionSign();
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

    SmallLocalIndexType directionIndex() const
    { return directionIdx_; }

    std::int_least8_t directionSign() const
    { return outerNormalSign_; }

    std::size_t scvfIdxWithCommonEntity() const
    { return scvfIdxWithCommonEntity_; }

private:
    GlobalPosition center_;
    GlobalPosition ipGlobal_;
    std::array<SmallLocalIndexType, 2> localScvIndices_;
    std::array<GridIndexType, 2> globalScvIndices_;
    SmallLocalIndexType localScvfIdx_;
    Scalar area_;
    SmallLocalIndexType directionIdx_;
    std::int_least8_t outerNormalSign_;
    GridIndexType globalScvfIdx_;
    GridIndexType scvfIdxWithCommonEntity_;
    FaceType faceType_;
    bool boundary_;
};


} // end namespace Dumux

#endif // DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUMEFACE_HH
