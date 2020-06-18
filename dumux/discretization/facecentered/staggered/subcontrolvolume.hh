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

#include <dumux/common/indextraits.hh>


#include <typeinfo>

namespace Dumux {

template <class GridView>
class FaceCenteredStaggeredSubControlVolume
{
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition =typename Element::Geometry::GlobalCoordinate;
    using Scalar = typename GridView::ctype;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SmallLocalIndexType = typename IndexTraits<GridView>::SmallLocalIndex;

public:

    FaceCenteredStaggeredSubControlVolume() = default;

    FaceCenteredStaggeredSubControlVolume(const GlobalPosition& center,
                                          const GlobalPosition& dofPosition,
                                          const Scalar volume,
                                          const GridIndexType globalIndex,
                                          const SmallLocalIndexType indexInElement,
                                          const GridIndexType dofIdx,
                                          const SmallLocalIndexType dirIdx,
                                          const std::int_least8_t dirSign,
                                          const GridIndexType eIdx,
                                          const bool boundary)
    : center_(center)
    , dofPosition_(dofPosition)
    , volume_(volume)
    , globalIndex_(globalIndex)
    , indexInElement_(indexInElement)
    , dofIdx_(dofIdx)
    , directionIdx_(dirIdx)
    , directionSign_(dirSign)
    , eIdx_(eIdx)
    , boundary_(boundary) {}

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
    { return eIdx_; }

    SmallLocalIndexType indexInElement() const
    { return indexInElement_; }

    SmallLocalIndexType localDofIndex() const
    { return indexInElement_; }

    SmallLocalIndexType directionIndex() const
    { return directionIdx_; }

    std::int_least8_t directionSign() const
    { return directionSign_; }

    bool boundary() const
    { return boundary_; }

private:
    GlobalPosition center_;
    GlobalPosition dofPosition_;
    Scalar volume_;
    GridIndexType globalIndex_;
    SmallLocalIndexType indexInElement_;
    GridIndexType dofIdx_;
    SmallLocalIndexType directionIdx_;
    std::int_least8_t directionSign_;
    GridIndexType eIdx_;
    bool boundary_;
};


} // end namespace Dumux

#endif // DUMUX_DISCRETIZATION_FACECENTERED_STAGGERED_SUBCONTROLVOLUME_HH
