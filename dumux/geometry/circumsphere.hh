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
 * \brief A function to compute the circumsphere of a given geometry
 * \note The circumscribed sphere (or circumsphere) is a sphere touching all vertices of the geometry
 *       Therefore it does not necessarily exist for all geometries.
 */
#ifndef DUMUX_GEOMETRY_CIRCUMSPHERE_HH
#define DUMUX_GEOMETRY_CIRCUMSPHERE_HH

#include <type_traits>
#include <dumux/common/math.hh>
#include <dumux/geometry/sphere.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Computes the circumsphere of a triangle
 * \note See https://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
 */
template<class Point>
static inline Sphere<typename Point::value_type, Point::dimension>
circumSphereOfTriangle(const Point& a, const Point& b, const Point& c)
{
    const auto ac = c - a;
    const auto ab = b - a;
    const auto n = crossProduct(ab, ac);
    const auto distCenterToA = (crossProduct(n, ab)*ac.two_norm2() + crossProduct(ac, n)*ab.two_norm2()) / (2.0*n.two_norm2());

    return { a + distCenterToA, distCenterToA.two_norm() };
}

/*!
 * \ingroup Geometry
 * \brief Computes the circumsphere of a triangle
 */
template<class Geometry, typename std::enable_if_t<Geometry::mydimension == 2, int> = 0>
static inline Sphere<typename Geometry::ctype, Geometry::coorddimension>
circumSphereOfTriangle(const Geometry& triangle)
{
    assert(triangle.corners() == 3 && "Geometry is not a triangle.");
    return circumSphereOfTriangle(triangle.corner(0), triangle.corner(1), triangle.corner(2));
}

} // end namespace Dumux

#endif
