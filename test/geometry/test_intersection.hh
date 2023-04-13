//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Common
 * \brief Intersection test helpers
 */
#ifndef DUMUX_TEST_INTERSECTION_HH
#define DUMUX_TEST_INTERSECTION_HH

#include <dumux/geometry/intersectspointgeometry.hh>

namespace Dumux {

/*!
 * \file
 * \ingroup Common
 * \brief Test intersection overload for point and geometry
 */
template<class Geometry>
bool testIntersection(const Geometry& geo,
                      const typename Geometry::GlobalCoordinate& p,
                      bool foundExpected, bool verbose)
{
    const bool found = intersectsPointGeometry(p, geo);
    if (!found && foundExpected)
    {
        std::cerr << "  Failed detecting intersection of " << geo.type();
        for (int i = 0; i < geo.corners(); ++i)
            std::cerr << " (" << geo.corner(i) << ")";
        std::cerr << " with point: " << p << std::endl;
    }
    else if (found && !foundExpected)
    {
        std::cerr << "  Found false positive: intersection of " << geo.type();
        for (int i = 0; i < geo.corners(); ++i)
            std::cerr << " (" << geo.corner(i) << ")";
        std::cerr << " with point: " << p << std::endl;
    }
    if (verbose)
    {
        if (found && foundExpected)
            std::cout << "  Found intersection with " << p << std::endl;
        else if (!found && !foundExpected)
            std::cout << "  No intersection with " << p << std::endl;
    }
    return (found == foundExpected);
}

} // end namespace Dumux

#endif
