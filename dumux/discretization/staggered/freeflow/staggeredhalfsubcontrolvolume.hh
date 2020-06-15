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
 * \copydoc Dumux::FreeFlowStaggeredHalfSubControlVolume
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_STAGGERED_HALF_SUBCONTROLVOLUME_HH
#define DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_STAGGERED_HALF_SUBCONTROLVOLUME_HH

#include <array>
#include <utility>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

#include <dumux/common/indextraits.hh>
#include <dumux/common/typetraits/isvalid.hh>


#include <typeinfo>

namespace Dumux {

template <class GridView>
class FreeFlowStaggeredHalfSubControlVolume
{
    using Grid = typename GridView::Grid;
    static constexpr int dim = Grid::dimension;
    static constexpr int dimWorld = Grid::dimensionworld;
    using Scalar = double;  //TODO

    using GlobalPosition = Dune::FieldVector<typename GridView::ctype, dim>;
public:

    FreeFlowStaggeredHalfSubControlVolume() = default;

    FreeFlowStaggeredHalfSubControlVolume(const GlobalPosition& center,
                                          const GlobalPosition& dofPosition,
                                          const Scalar volume,
                                          const std::size_t dofIdx,
                                          const std::size_t correspondingCellCenterScvfIdx, // TODO remove or rename
                                          const std::size_t dirIdx,
                                          const Scalar dirSign,
                                          const bool boundary)
    : center_(center)
    , dofPosition_(dofPosition)
    , volume_(volume)
    , dofIdx_(dofIdx)
    , correspondingCellCenterScvfIdx_(correspondingCellCenterScvfIdx)
    , directionIdx_(dirIdx)
    , directionSign_(dirSign)
    , boundary_(boundary) {}

    //! The center of the sub control volume
    const GlobalPosition& center() const
    { return center_; }

    //! The position of the degree of freedom
    const GlobalPosition& dofPosition() const
    { return dofPosition_; }

    Scalar volume() const
    { return volume_; }

    std::size_t dofIndex() const
    { return dofIdx_; }

    std::size_t correspondingCellCenterScvfIndex() const
    { return correspondingCellCenterScvfIdx_; }

    std::size_t directionIndex() const
    { return directionIdx_; }

    Scalar directionSign() const
    { return directionSign_; }

    bool boundary() const
    { return boundary_; }

private:
    GlobalPosition center_;
    GlobalPosition dofPosition_;
    Scalar volume_;
    std::size_t dofIdx_;
    std::size_t correspondingCellCenterScvfIdx_;
    std::size_t directionIdx_;
    Scalar directionSign_;
    bool boundary_;
};


} // end namespace Dumux

#endif // DUMUX_DISCRETIZATION_STAGGERED_FREE_FLOW_SUBCONTROLVOLUMEFACE_HH
