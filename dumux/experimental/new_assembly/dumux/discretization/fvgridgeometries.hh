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
 * \ingroup FVDiscretization
 * \brief Default implementations of the grid geometries for finite-volume schemes.
 */
#ifndef DUMUX_DISCRETIZATION_FV_GRID_GEOMETRIES_HH
#define DUMUX_DISCRETIZATION_FV_GRID_GEOMETRIES_HH

#include <concepts>
#include <utility>
#include <type_traits>

namespace Dumux {

/*!
 * \ingroup FVDiscretization
 * \brief Represents a face between two control volumes.
 */
template<typename Coordinate,
         typename ctype = typename Coordinate::value_type>
class Face
{
public:
    Face() = default;

    template<typename C> requires(std::same_as<Coordinate, std::decay_t<C>>)
    Face(ctype area, C&& center)
    : area_(area)
    , center_(std::forward<C>(center))
    {}

    ctype area() const { return area_; }
    const Coordinate& center() const { return center_; }

private:
    ctype area_;
    Coordinate center_;
};

/*!
 * \ingroup FVDiscretization
 * \brief A sub-control volume face is a face with a defined orientation of the outer normal vector.
 *        That is, there exists a unique sub-control volume face for each neighboring control volume
 *        of a face.
 * \note The lifetime of this object is bound to the lifetime of the underlying `Face`
 */
template<typename Coordinate,
         typename ctype = typename Coordinate::value_type,
         typename Face = Dumux::Face<Coordinate, ctype>>
class SubControlVolumeFace
{
public:
    SubControlVolumeFace() = default;

    template<typename C> requires(std::same_as<Coordinate, std::decay_t<C>>)
    SubControlVolumeFace(const Face& face, C&& normal)
    : face_(&face)
    , normal_(std::forward<C>(normal))
    {}

    ctype area() const { return face_->area(); }
    const Coordinate& center() const { return face_->center(); }
    const Coordinate& ipGlobal() const { return face_->center(); }
    const Coordinate& unitOuterNormal() const { return normal_; }

private:
    const Face* face_{nullptr};
    Coordinate normal_;
};

/*!
 * \ingroup FVDiscretization
 * \brief Sub-control volume implementation for finite-volume schemes.
 */
template<typename GridIndex,
         typename Coordinate,
         typename ctype = typename Coordinate::value_type>
class SubControlVolume
{
public:
    SubControlVolume() = default;

    template<typename C> requires(std::same_as<Coordinate, std::decay_t<C>>)
    SubControlVolume(GridIndex dofIndex, ctype volume, C&& center)
    : dofIndex_(dofIndex)
    , volume_(volume)
    , center_(std::forward<C>(center))
    {}

    ctype volume() const { return volume_; }
    const Coordinate& center() const { return center_; }
    const Coordinate& dofPosition() const { return center_; }
    GridIndex dofIndex() const { return dofIndex_; }

private:
    GridIndex dofIndex_;
    ctype volume_;
    Coordinate center_;
};

} // end namespace Dumux

#endif
