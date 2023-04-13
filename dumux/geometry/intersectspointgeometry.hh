// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief Detect if a point intersects a geometry
 */
#ifndef DUMUX_GEOMETRY_INTERSECTS_POINT_GEOMETRY_HH
#define DUMUX_GEOMETRY_INTERSECTS_POINT_GEOMETRY_HH

#include <cmath>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dumux/geometry/intersectspointsimplex.hh>

namespace Dumux {

/*!
 * \ingroup Geometry
 * \brief Find out whether a point is inside a three-dimensional geometry
 */
template <class ctype, int dimworld, class Geometry, typename std::enable_if_t<(Geometry::mydimension == 3), int> = 0>
bool intersectsPointGeometry(const Dune::FieldVector<ctype, dimworld>& point, const Geometry& g)
{
    // get the g type
    const auto type = g.type();

    // if it's a tetrahedron we can check directly
    if (type.isTetrahedron())
        return intersectsPointSimplex(point, g.corner(0), g.corner(1), g.corner(2), g.corner(3));

    // split hexahedrons into five tetrahedrons
    else if (type.isHexahedron())
    {
        if (intersectsPointSimplex(point, g.corner(0), g.corner(1), g.corner(3), g.corner(5))) return true;
        if (intersectsPointSimplex(point, g.corner(0), g.corner(5), g.corner(6), g.corner(4))) return true;
        if (intersectsPointSimplex(point, g.corner(5), g.corner(3), g.corner(6), g.corner(7))) return true;
        if (intersectsPointSimplex(point, g.corner(0), g.corner(3), g.corner(2), g.corner(6))) return true;
        if (intersectsPointSimplex(point, g.corner(5), g.corner(3), g.corner(0), g.corner(6))) return true;
        return false;
    }

    // split pyramids into two tetrahedrons
    else if (type.isPyramid())
    {
        if (intersectsPointSimplex(point, g.corner(0), g.corner(1), g.corner(2), g.corner(4))) return true;
        if (intersectsPointSimplex(point, g.corner(1), g.corner(3), g.corner(2), g.corner(4))) return true;
        return false;
    }

    // split prisms into three tetrahedrons
    else if (type.isPrism())
    {
        if (intersectsPointSimplex(point, g.corner(0), g.corner(1), g.corner(2), g.corner(4))) return true;
        if (intersectsPointSimplex(point, g.corner(3), g.corner(0), g.corner(2), g.corner(4))) return true;
        if (intersectsPointSimplex(point, g.corner(2), g.corner(5), g.corner(3), g.corner(4))) return true;
        return false;
    }

    else
        DUNE_THROW(Dune::NotImplemented,
                   "Intersection for point and geometry type "
                   << type << " in " << dimworld << "-dimensional world.");
}

/*!
 * \ingroup Geometry
 * \brief Find out whether a point is inside a two-dimensional geometry
 */
template <class ctype, int dimworld, class Geometry, typename std::enable_if_t<(Geometry::mydimension == 2), int> = 0>
bool intersectsPointGeometry(const Dune::FieldVector<ctype, dimworld>& point, const Geometry& g)
{
    // get the g type
    const auto type = g.type();

    // if it's a triangle we can check directly
    if (type.isTriangle())
        return intersectsPointSimplex(point, g.corner(0), g.corner(1), g.corner(2));

    // split quadrilaterals into two triangles
    else if (type.isQuadrilateral())
    {
        if (intersectsPointSimplex(point, g.corner(0), g.corner(1), g.corner(3))) return true;
        if (intersectsPointSimplex(point, g.corner(0), g.corner(3), g.corner(2))) return true;
        return false;
    }

    else
        DUNE_THROW(Dune::NotImplemented,
                   "Intersection for point and geometry type "
                   << type << " in " << dimworld << "-dimensional world.");
}

/*!
 * \ingroup Geometry
 * \brief Find out whether a point is inside a one-dimensional geometry
 */
template <class ctype, int dimworld, class Geometry, typename std::enable_if_t<(Geometry::mydimension == 1), int> = 0>
bool intersectsPointGeometry(const Dune::FieldVector<ctype, dimworld>& point, const Geometry& g)
{
    return intersectsPointSimplex(point, g.corner(0), g.corner(1));
}

} // end namespace Dumux

#endif
