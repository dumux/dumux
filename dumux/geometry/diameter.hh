// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
