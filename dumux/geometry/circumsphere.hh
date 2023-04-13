// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
