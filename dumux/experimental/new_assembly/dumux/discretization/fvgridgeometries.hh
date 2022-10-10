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

#include <utility>

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
    Face(ctype area, Coordinate&& center)
    : area_(area)
    , center_(std::move(center))
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
 */
template<typename Coordinate,
         typename ctype = typename Coordinate::value_type,
         typename Face = Dumux::Face<Coordinate, ctype>>
class SubControlVolumeFace
{
public:
    SubControlVolumeFace() = default;
    SubControlVolumeFace(const Face& face, Coordinate&& normal)
    : face_(&face)
    , normal_(std::move(normal))
    {}

    /*!
     * \brief Forbid copies to avoid issues with the lifetime of this
              object being bound to that of the underlying face
     */
    SubControlVolumeFace(const SubControlVolumeFace&) = delete;

    ctype area() const { return face_->area(); }
    const Coordinate& center() const { return face_->center(); }
    const Coordinate& ipGlobal() const { return face_->center(); }
    const Coordinate& unitOuterNormal() const { return normal_; }

private:
    const Face* face_{nullptr};
    Coordinate normal_;
};

//! Base class for sub-control volumes, storing basic geometric information.
template<typename GridIndex,
         typename Coordinate,
         typename ctype = typename Coordinate::value_type>
class SubControlVolumeBase
{
public:
    SubControlVolumeBase() = default;
    SubControlVolumeBase(GridIndex dofIndex,
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

private:
    GridIndex dofIndex_;
    ctype volume_;
    Coordinate center_;
};

/*!
 * \ingroup FVDiscretization
 * \brief Sub-control volume implementation for finite-volume schemes.
 */
template<typename GridIndex,
         typename Coordinate,
         typename ctype = typename Coordinate::value_type>
class SubControlVolume
: public SubControlVolumeBase<GridIndex, Coordinate, ctype>
{
    using ParentType = SubControlVolumeBase<GridIndex, Coordinate, ctype>;

public:
    SubControlVolume() = default;
    SubControlVolume(GridIndex dofIndex,
                     ctype volume,
                     Coordinate&& center,
                     unsigned int localIndex)
    : ParentType(dofIndex, volume, std::move(center))
    , localIndex_(localIndex)
    {}

    unsigned int localDofIndex() const { return localIndex_; }

private:
    unsigned int localIndex_;
};

/*!
 * \ingroup FVDiscretization
 * \brief Sub-control volume implementation for cell-centered finite-volume schemes.
 *        In cell-centered schemes, there is only a single sub-control volume per grid
 *        cell, and thus, the `localDofIndex` is zero and is known at compile-time.
 */
template<typename GridIndex,
         typename Coordinate,
         typename ctype = typename Coordinate::value_type>
class CCSubControlVolume : public SubControlVolumeBase<GridIndex, Coordinate, ctype>
{
    using ParentType = SubControlVolumeBase<GridIndex, Coordinate, ctype>;

public:
    using ParentType::ParentType;
    unsigned int localDofIndex() const { return 0; }
};

} // end namespace Dumux

#endif
