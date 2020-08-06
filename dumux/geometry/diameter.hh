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
 * \brief A function to compute a geometry's diameter, i.e.
 *        the longest distance between points of a geometry
 */
#ifndef DUMUX_GEOMETRY_DIAMETER_HH
#define DUMUX_GEOMETRY_DIAMETER_HH

#include <algorithm>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Computes the longest distance between points of a geometry
 * \note Useful e.g. to compute the maximum cell diameter of a grid
 */
template<class Geometry>
typename Geometry::ctype diameter(const Geometry& geo)
{
    using std::max;
    typename Geometry::ctype h = 0.0;
    for (std::size_t i = 0; i < geo.corners(); ++i)
        for (std::size_t j = i + 1; j < geo.corners(); ++j)
            h = max(h, (geo.corner(i)-geo.corner(j)).two_norm());

    return h;
}

} // end namespace Dumux

#endif
