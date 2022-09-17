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
 * \ingroup CCTpfaDiscretization
 * \brief Default implementations of the grid geometries for cell-centered schemes.
 */
#ifndef DUMUX_DISCRETIZATION_CC_GRID_GEOMETRIES_HH
#define DUMUX_DISCRETIZATION_CC_GRID_GEOMETRIES_HH

#include <utility>

namespace Dumux {

template<typename Coordinate,
         typename ctype = typename Coordinate::value_type>
class CCFace
{
public:
    CCFace() = default;
    CCFace(ctype area, Coordinate&& center)
    : area_(area)
    , center_(std::move(center))
    {}

    ctype area() const { return area_; }
    const Coordinate& center() const { return center_; }

private:
    ctype area_;
    Coordinate center_;
};

template<typename Coordinate,
         typename ctype = typename Coordinate::value_type,
         typename Face = CCFace<Coordinate, ctype>>
class CCSubControlVolumeFace
{
public:
    CCSubControlVolumeFace() = default;
    CCSubControlVolumeFace(const Face& face, Coordinate&& normal)
    : face_(&face)
    , normal_(std::move(normal))
    {}

    ctype area() const { return face_->area(); }
    const Coordinate& center() const { return face_->center(); }
    const Coordinate& ipGlobal() const { return face_->center(); }
    const Coordinate& unitOuterNormal() const { return normal_; }

private:
    const Face* face_{nullptr};
    Coordinate normal_;
};

template<typename GridIndex,
         typename Coordinate,
         typename ctype = typename Coordinate::value_type>
class CCSubControlVolume
{
public:
    CCSubControlVolume() = default;
    CCSubControlVolume(GridIndex dofIndex,
                       ctype volume,
                       Coordinate&& center)
    : dofIndex_(dofIndex)
    , volume_(volume)
    , center_(std::move(center))
    {}

    ctype volume() const { return volume_; }
    const Coordinate& center() const { return center_; }
    const Coordinate& dofPosition() const { return center_; }

    GridIndex dofIndex() const { return dofIndex_; }
    unsigned int localDofIndex() const { return 0; }

private:
    GridIndex dofIndex_;
    ctype volume_;
    Coordinate center_;
};

} // end namespace Dumux

#endif
