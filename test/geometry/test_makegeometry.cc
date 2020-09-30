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
 *
 * \brief Test for the makegeometry method
 */
 #include <config.h>

 #include <iostream>
 #include <algorithm>
 #include <random>

#include <dune/common/float_cmp.hh>
#include <dumux/geometry/intersectspointgeometry.hh>
#include <dumux/geometry/makegeometry.hh>
#include <dumux/geometry/grahamconvexhull.hh>

//! test if the ordering of the input points has an effect on the resulting geometry
template<class GlobalPosition>
void permutatePointsAndTest(const std::vector<GlobalPosition>& cornerPoints,
                            const std::vector<GlobalPosition>& pointsWithinGeometry,
                            const std::vector<GlobalPosition>& pointsOutsideGeometry,
                            const bool verbose = false)
{
    std::array<int, 4> s = {0,1,2,3};

    if (Dumux::grahamConvexHull<2>(cornerPoints).size() != 4)
    {
        if (verbose)
            std::cout << "Input points do not span a strictly convex polygon. Skipping." << std::endl;
        return;
    }

    const auto area = Dumux::makeDuneQuadrilaterial(cornerPoints).volume();

    do {
            const std::vector<GlobalPosition> corners = {cornerPoints[s[0]], cornerPoints[s[1]], cornerPoints[s[2]], cornerPoints[s[3]]};

            if (Dumux::grahamConvexHull<2>(corners).size() != 4)
            {
                if (verbose)
                    std::cout << "Input points do not span a strictly convex polygon. Skipping." << std::endl;
                continue;
            }

            const auto quad = Dumux::makeDuneQuadrilaterial(corners);

            auto printCorners = [&quad]()
            {
                std::ostringstream tmp;
                for (int i = 0; i < quad.corners(); ++i)
                    tmp << "(" << quad.corner(i) << ")  ";
                return tmp.str();
            };

            if (verbose)
            {
                std::cout << "resulting corners \n";
                std::cout << printCorners();
                std::cout << std::endl;

                std::cout << "actual area: " << area << ", area after permuation: " << quad.volume() << std::endl;
            }

            if (!Dune::FloatCmp::eq(quad.volume(), area))
                DUNE_THROW(Dune::InvalidStateException, "Area of quadrilateral after permuation of input points is wrong");

            if (verbose)
                std::cout  << "checking point intersection: \n";

            for (const auto& p: pointsWithinGeometry)
            {
                if (Dumux::intersectsPointGeometry(p, quad))
                {
                    if(verbose)
                        std::cout << "point " << p << " lies within the quadrilateral" << std::endl;
                }
                else
                    DUNE_THROW(Dune::InvalidStateException, "False negative: Check for point " << p << " which is inside the geometry " << printCorners() << " failed");
            }

            for (const auto& p : pointsOutsideGeometry)
            {
                if (!Dumux::intersectsPointGeometry(p, quad))
                {
                    if (verbose)
                        std::cout << "point " << p << " lies outside of the quadrilateral" << std::endl;
                }
                else
                    DUNE_THROW(Dune::InvalidStateException, "False positive: Check for point " << p << " which is outside the geometry " << printCorners() << "  failed");
            }

    } while (std::next_permutation(s.begin(), s.end()));
}

template<class GlobalPosition>
void checkAxisAlignedGeometry(std::vector<GlobalPosition>& cornerPoints,
                              std::vector<GlobalPosition>& pointsWithinGeometry,
                              std::vector<GlobalPosition>& pointsOutsideGeometry,
                              const int normalDirection,
                              const typename GlobalPosition::value_type scale)
{
    static const char dim[]= "xyz";
    std::cout << "testing for quadrilateral with normal in " << dim[normalDirection] << " direction" << std::endl;

    // check if points are coplanar
    if(!Dumux::pointsAreCoplanar(cornerPoints))
        DUNE_THROW(Dune::InvalidStateException, "False negative, points are actually coplanar!");

    cornerPoints[0][normalDirection] += 1e-5*scale; // we make them non-coplanar
    if(Dumux::pointsAreCoplanar(cornerPoints))
        DUNE_THROW(Dune::InvalidStateException, "False positive, points are actually not coplanar!");

    cornerPoints[0][normalDirection] -= 1e-5*scale;

    permutatePointsAndTest(cornerPoints, pointsWithinGeometry, pointsOutsideGeometry);

    const auto origCornerPoints = cornerPoints;

    std::vector<GlobalPosition> pointsWithinGeometryForRandomTest = { {0.6, 0.6, 0.6},
                                                                     {0.4, 0.4, 0.4} };
    for(auto& p : pointsWithinGeometryForRandomTest)
    {
        p *= scale;
        p[normalDirection] = 0.0;
    }

    auto pointsOutsideGeometryForRandomTest = pointsOutsideGeometry;
    pointsOutsideGeometryForRandomTest[0] = {1.5*scale, 1.5*scale, 1.5*scale};
    pointsOutsideGeometryForRandomTest[0][normalDirection] = 0.0;

    // uniform random number generator
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<> uniformdist(-0.3*scale, 0.3*scale);

    for(int i = 0; i < 10; i++)
    {
        for(auto&& p : cornerPoints)
        {
            for(int x = 0; x < p.size(); ++x)
            {
                if(x != normalDirection)
                    p[x] += uniformdist(generator);
            }
        }

        permutatePointsAndTest(cornerPoints, pointsWithinGeometryForRandomTest, pointsOutsideGeometryForRandomTest);

        cornerPoints = origCornerPoints;
    }
}

int main(int argc, char** argv)
{
    using namespace Dumux;
    using GlobalPosition = Dune::FieldVector<double, 3>;

    std::array<double, 3> scaling{{1e-12, 1.0, 1e12}};

    for (const double scale : scaling)
    {
        const double size = 1.0*scale;
        const double half = 0.5*scale;
        const double small = 1e-3*scale;

        GlobalPosition p0 = {0, 0, 0};
        GlobalPosition p1 = {size, 0, 0};
        GlobalPosition p2 = {0, size, 0};
        GlobalPosition p3 = {size, size, 0};

        std::vector<GlobalPosition> cornerPoints = {p0, p1, p2, p3};

        std::vector<GlobalPosition> pointsWithinGeometry = { GlobalPosition{half, half, 0.0},
                                                             GlobalPosition{cornerPoints[0][0] + small, cornerPoints[0][1] + small, 0.0},
                                                             GlobalPosition{cornerPoints[3][0] - small, cornerPoints[3][1] - small, 0.0} };

        std::vector<GlobalPosition> pointsOutsideGeometry = { GlobalPosition{cornerPoints[0][0] - small, cornerPoints[0][1] - small, 0.0},
                                                              GlobalPosition{half, half, small} };

        // do the checks for a quadrilateral parallel to the x and y axis
        checkAxisAlignedGeometry(cornerPoints, pointsWithinGeometry, pointsOutsideGeometry, 2, size);

        // rotate the quadrilateral to make it parallel to the other axes and test again
        for(int i = 1; i >=0; --i)
        {
            for(auto& p : cornerPoints)
                std::rotate(p.begin(), p.begin() + 1, p.end());
            for(auto& p : pointsWithinGeometry)
                std::rotate(p.begin(), p.begin() + 1, p.end());
            for(auto& p : pointsOutsideGeometry)
                std::rotate(p.begin(), p.begin() + 1, p.end());

            checkAxisAlignedGeometry(cornerPoints, pointsWithinGeometry, pointsOutsideGeometry, i, size);
        }

        std::cout << "testing for non axis-aligned quadrilateral" << std::endl;

        cornerPoints[0][0] += half;
        cornerPoints[1][0] -= half;
        cornerPoints[2][0] += half;
        cornerPoints[3][0] -= half;

        GlobalPosition pointToCheck5 = {0.0, half, half};

        pointsWithinGeometry = {pointToCheck5};

        pointsOutsideGeometry[1] = pointsOutsideGeometry[0];

        permutatePointsAndTest(cornerPoints, pointsWithinGeometry, pointsOutsideGeometry);
    }

    return 0;
}
