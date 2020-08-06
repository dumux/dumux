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
