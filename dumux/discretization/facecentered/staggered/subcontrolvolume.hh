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
 * \copydoc Dumux::FaceCenteredStaggeredSubControlVolume
 */
#ifndef DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUME_HH

#include <array>
#include <utility>
#include <dune/geometry/type.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>

#include <dumux/common/indextraits.hh>


#include <typeinfo>

namespace Dumux {

/*!
 * \ingroup CCDiscretization
 * \brief Default traits class to be used for the sub-control volumes
 *        for the cell-centered finite volume scheme using TPFA
 * \tparam GV the type of the grid view
 */
template<class GridView>
struct FaceCenteredDefaultScvGeometryTraits
{
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using LocalIndexType = typename IndexTraits<GridView>::LocalIndex;
    using Scalar = typename GridView::ctype;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dim = GridView::Grid::dimension;
    static constexpr int dimWorld = GridView::Grid::dimensionworld;
    using CornerStorage = std::array<GlobalPosition, (1<<(dim))>;
    using Geometry = Dune::AxisAlignedCubeGeometry<Scalar, dim, dimWorld>;
};


template<class GridView, class T = FaceCenteredDefaultScvGeometryTraits<GridView>>
class FaceCenteredStaggeredSubControlVolume
{
    using Geometry = typename T::Geometry;
    using CornerStorage = typename T::CornerStorage;
    using Element = typename T::Element;
    using GlobalPosition = typename T::GlobalPosition;
    using Scalar = typename T::Scalar;
    using GridIndexType = typename T::GridIndexType;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;

public:
    //! state the traits public and thus export all types
    using Traits = T;

    FaceCenteredStaggeredSubControlVolume() = default;

    template<class Corners>
    FaceCenteredStaggeredSubControlVolume(const GlobalPosition& center,
                                          const GlobalPosition& dofPosition,
                                          Corners&& corners,
                                          const Scalar volume,
                                          const GridIndexType globalIndex,
                                          const SmallLocalIndexType indexInElement,
                                          const GridIndexType dofIdx,
                                          const SmallLocalIndexType dofAxis,
                                          const std::int_least8_t dirSign,
                                          const GridIndexType elementIdx,
                                          const GridIndexType neighborElementIdx,
                                          const bool boundary)
    : center_(center)
    , dofPosition_(dofPosition)
    , volume_(volume)
    , globalIndex_(globalIndex)
    , indexInElement_(indexInElement)
    , dofIdx_(dofIdx)
    , dofAxis_(dofAxis)
    , directionSign_(dirSign)
    , elementIdx_(elementIdx)
    , neighborElementIdx_(neighborElementIdx)
    , boundary_(boundary)
    {
        if constexpr (std::is_same_v<Corners, CornerStorage>)
            corners_ = std::move(corners);
        else
        {
            for (int i = 0; i < corners_.size(); ++i)
                corners_[i] = corners[i];
        }
    }

    //! The center of the sub control volume
    const GlobalPosition& center() const
    { return center_; }

    //! The position of the degree of freedom
    const GlobalPosition& dofPosition() const
    { return dofPosition_; }

    Scalar volume() const
    { return volume_; }

    GridIndexType dofIndex() const
    { return dofIdx_; }

    GridIndexType index() const
    { return globalIndex_; }

    GridIndexType elementIndex() const
    { return elementIdx_; }

    GridIndexType neighborElementIdx() const
    { return neighborElementIdx_; }

    SmallLocalIndexType indexInElement() const
    { return indexInElement_; }

    SmallLocalIndexType localDofIndex() const
    { return indexInElement_; }

    SmallLocalIndexType dofAxis() const
    { return dofAxis_; }

    std::int_least8_t directionSign() const
    { return directionSign_; }

    bool boundary() const
    { return boundary_; }

    const GlobalPosition& corner(unsigned int localIdx) const
    {
        assert(localIdx < corners_.size() && "provided index exceeds the number of corners");
        return corners_[localIdx];
    }

    //! The geometry of the sub control volume face
    Geometry geometry() const
    {
        return Geometry(corners_.front(), corners_.back());
    }

private:
    GlobalPosition center_;
    GlobalPosition dofPosition_;
    CornerStorage corners_;
    Scalar volume_;
    GridIndexType globalIndex_;
    SmallLocalIndexType indexInElement_;
    GridIndexType dofIdx_;
    SmallLocalIndexType dofAxis_;
    std::int_least8_t directionSign_;
    GridIndexType elementIdx_;
    GridIndexType neighborElementIdx_;
    bool boundary_;
};

} // end namespace Dumux

#endif // DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUME_HH
