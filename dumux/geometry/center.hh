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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Geometry
 * \brief Compute the center point of a convex polytope geometry or a random-access container of corner points
 */
#ifndef DUMUX_GEOMETRY_CENTER_HH
#define DUMUX_GEOMETRY_CENTER_HH

#include <iterator>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief The center of a given list of corners
 */
template<class Corners>
typename Corners::value_type center(const Corners& corners)
{
    auto it = std::begin(corners);
    auto center = *it; ++it;
    for (; it != std::end(corners); ++it)
        center += *it;
    center /= std::distance(std::begin(corners), std::end(corners));
    return center;
}

} // end namespace Dumux

#endif
