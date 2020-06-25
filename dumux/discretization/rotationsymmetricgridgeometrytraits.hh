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
 * \ingroup Discretization
 * \brief Rotation policy for defining rotational symmetric grid geometries
 */
#ifndef DUMUX_DISCRETIZATION_ROTATION_SYMMETRIC_GG_TRAITS_HH
#define DUMUX_DISCRETIZATION_ROTATION_SYMMETRIC_GG_TRAITS_HH

#warning "This header is deprecated and will be removed after release 3.3"
#include <dumux/discretization/rotationpolicy.hh>
#include <dumux/discretization/rotationsymmetricscv.hh>
#include <dumux/discretization/rotationsymmetricscvf.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Traits for rotation symmetric grid geometries
 * \tparam BaseTraits The traits type to turn into rotation symmetric traits
 * \tparam rotPolicy The rotation policy (see RotationPolicy enum class)
 */
template<class BaseTraits, RotationPolicy rotPolicy>
struct [[deprecated("Will be removed after release 3.3. Use Extrusion from extrusion.hh")]] RotationSymmetricGridGeometryTraits : public BaseTraits
{
    using SubControlVolume = RotationSymmetricSubControlVolume<typename BaseTraits::SubControlVolume, rotPolicy>;
    using SubControlVolumeFace = RotationSymmetricSubControlVolumeFace<typename BaseTraits::SubControlVolumeFace, rotPolicy>;
};

} // end namespace Dumux

#endif
