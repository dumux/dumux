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
#ifndef DUMUX_DISCRETIZATION_ROTATION_POLICY_HH
#define DUMUX_DISCRETIZATION_ROTATION_POLICY_HH

#warning "This header is deprecated and will be removed after release 3.3"

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Rotation policies
 * - disc (or annulus): rotate a segment around an axis perpendicular to the line through the segment (we choose )
 * - ball (or shell): rotate a segment around a point on the line through the segment
 * - toroid: rotate a polygon around one of its axis (we rotate about the second axis)
 */
enum class [[deprecated("Will be removed after release 3.3. Use Extrusion from extrusion.hh")]] RotationPolicy
{ disc, ball, toroid };

} // end namespace Dumux

#endif
