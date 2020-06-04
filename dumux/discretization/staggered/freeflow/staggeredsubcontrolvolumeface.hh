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
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_STAGGEREDSUBCONTROLVOLUMEFACE_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_STAGGEREDSUBCONTROLVOLUMEFACE_HH

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
class FreeFlowRealStaggeredSubControlVolumeFace
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;

    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;

    using Scalar = double; // TODO traits


    using GlobalPosition = Dune::FieldVector<typename GridView::ctype, dim>;
public:
    enum class FaceType {frontal, lateral};

    FreeFlowRealStaggeredSubControlVolumeFace() = default;

    FreeFlowRealStaggeredSubControlVolumeFace(const GlobalPosition& center,
                                              const GlobalPosition& ipGlobal,
                                              const std::array<GridIndexType, 2> scvIndices,
                                              const Scalar area,
                                              const LocalIndexType directionIdx,
                                              const int outerNormalSign,
                                              const std::size_t localScvfIdx,
                                              const FaceType faceType,
                                              const bool boundary)
    : center_(center)
    , ipGlobal_(ipGlobal)
    , scvIndices_(scvIndices)
    , area_(area)
    , directionIdx_(directionIdx)
    , outerNormalSign_(outerNormalSign)
    , localScvfIdx_(localScvfIdx)
    , faceType_(faceType)
    , boundary_(boundary) {}

    //! The center of the sub control volume face
    const GlobalPosition& center() const
    { return center_; }

    //! The integration point of the sub control volume face
    const GlobalPosition& ipGlobal() const // TODO
    { return center_; }

    //! Index of the inside sub control volume for spatial param evaluation
    GridIndexType insideScvIdx() const
    { return scvIndices_[0]; }

    //! index of the outside sub control volume for spatial param evaluation
    GridIndexType outsideScvIdx() const
    { return scvIndices_[1]; }

    std::size_t localScvfIdx() const
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

    LocalIndexType directionIndex() const
    { return directionIdx_; }

    int directionSign() const // TODO smaller type
    { return outerNormalSign_; }

private:
    GlobalPosition center_;
    GlobalPosition ipGlobal_;
    std::array<GridIndexType, 2> scvIndices_; // TODO higher order?
    Scalar area_;
    LocalIndexType directionIdx_;
    int outerNormalSign_;
    std::size_t localScvfIdx_;
    FaceType faceType_;
    bool boundary_;
};


} // end namespace Dumux

#endif // DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_SUBCONTROLVOLUMEFACE_HH
