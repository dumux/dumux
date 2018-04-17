// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
#include <dumux/common/geometry/intersectspointgeometry.hh>
#include <dumux/common/geometry/makegeometry.hh>

//! test if the ordering of the input points has an effect on the resulting geometry
template<class GlobalPosition>
void permutatePointsAndTest(const std::vector<GlobalPosition>& cornerPoints,
                            const std::vector<GlobalPosition>& pointsWithinGeometry,
                            const std::vector<GlobalPosition>& pointsOutsideGeometry,
                            const bool verbose = false)
{
    std::array<int, 4> s = {0,1,2,3};

    const auto area = Dumux::makeDuneQuadrilaterial(cornerPoints).volume();

    do {
            std::vector<GlobalPosition> tmp = {cornerPoints[s[0]], cornerPoints[s[1]], cornerPoints[s[2]], cornerPoints[s[3]]};

            if(verbose)
            {
                std::cout << "input corner points: " << std::endl;
                for(auto&& i : tmp)
                    std::cout << "(" << i << ")  ";
                std::cout << std::endl;
            }

            const auto quad = Dumux::makeDuneQuadrilaterial(tmp);

            if(verbose)
            {
                std::cout << "resulting corners \n";
                for(int i = 0; i < quad.corners(); ++i)
                    std::cout << "(" << quad.corner(i) << ")  ";
                std::cout << std::endl;

                std::cout << "actual area: " << area << ", area after permuation: " << quad.volume() << std::endl;
            }

            if(!Dune::FloatCmp::eq(quad.volume(), area))
                DUNE_THROW(Dune::InvalidStateException, "Area of quadrilateral after permuation of input points is wrong");

            if(verbose)
                std::cout  << "checking point intersection: \n";

            for(const auto& p: pointsWithinGeometry)
            {
                if(Dumux::intersectsPointGeometry(p, quad))
                {
                    if(verbose)
                        std::cout << "point " << p << " lies within the quadrilateral" << std::endl;
                }
                else
                    DUNE_THROW(Dune::InvalidStateException, "Check for point inside geometry failed. Point " << p << " does not lie within the geometry!");
            }

            for(const auto& p : pointsOutsideGeometry)
            {
                if(!Dumux::intersectsPointGeometry(p, quad))
                {
                    if(verbose)
                        std::cout << "point " << p << " lies outside of the quadrilateral" << std::endl;
                }
                else
                    DUNE_THROW(Dune::InvalidStateException, "Check for point outside geometry failed. Point " << p << " does lie within the geometry!");
            }

    } while(std::next_permutation(s.begin(), s.end()));
}

template<class GlobalPosition>
void checkAxisAlignedGeometry(std::vector<GlobalPosition>& cornerPoints,
                              std::vector<GlobalPosition>& pointsWithinGeometry,
                              std::vector<GlobalPosition>& pointsOutsideGeometry,
                              const int normalDirection)
{
    static const char dim[]= "xyz";
    std::cout << "testing for quadrilateral with normal in " << dim[normalDirection] << " direction" << std::endl;

    // check if points are coplanar
    if(!Dumux::pointsAreCoplanar(cornerPoints))
        DUNE_THROW(Dune::InvalidStateException, "False positive, points are actually coplanar!");

    cornerPoints[0][normalDirection] += 1e-9;
    if(Dumux::pointsAreCoplanar(cornerPoints))
        DUNE_THROW(Dune::InvalidStateException, "Points are not coplanar!");

    cornerPoints[0][normalDirection] -= 1e-9;

    permutatePointsAndTest(cornerPoints, pointsWithinGeometry, pointsOutsideGeometry);

    const auto origCornerPoints = cornerPoints;

    std::vector<GlobalPosition> pointsWithinGeometryForRandomTest = { {0.6, 0.6, 0.6},
                                                                     {0.4, 0.4, 0.4} };
    for(auto& p : pointsWithinGeometryForRandomTest)
        p[normalDirection] = 0.0;

    auto pointsOutsideGeometryForRandomTest = pointsOutsideGeometry;
    pointsOutsideGeometryForRandomTest[0] = {1.5, 1.5, 1.5};
    pointsOutsideGeometryForRandomTest[0][normalDirection] = 0.0;

    for(int i = 0; i< 10; i++)
    {
        // uniform random number generator
        std::random_device rd;
        std::mt19937 generator(rd());
        std::uniform_real_distribution<> uniformdist(-0.3, 0.3);

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

int main(int argc, char** argv) try
{
    using namespace Dumux;
    using GlobalPosition = Dune::FieldVector<double, 3>;

    GlobalPosition p0 = {0,0,0};
    GlobalPosition p1 = {1,0,0};
    GlobalPosition p2 = {0,1,0};
    GlobalPosition p3 = {1,1,0};

    std::vector<GlobalPosition> cornerPoints = {p0, p1, p2, p3};

    std::vector<GlobalPosition> pointsWithinGeometry = { GlobalPosition{0.5, 0.5, 0.0},
                                                         GlobalPosition{cornerPoints[0][0] + 1e-3, cornerPoints[0][1] + 1e-3, 0.0},
                                                         GlobalPosition{cornerPoints[3][0] - 1e-3, cornerPoints[3][1] - 1e-3, 0.0} };

    std::vector<GlobalPosition> pointsOutsideGeometry = { GlobalPosition{cornerPoints[0][0] - 1e-3, cornerPoints[0][1] - 1e-3, 0.0},
                                                          GlobalPosition{0.5, 0.5, 1e-3} };

    // do the checks for a quadrilateral parallel to the x and y axis
    checkAxisAlignedGeometry(cornerPoints, pointsWithinGeometry, pointsOutsideGeometry, 2);

    // rotate the quadrilateral to make it parallel to the other axes and test again
    for(int i = 1; i >=0; --i)
    {
        for(auto& p : cornerPoints)
            std::rotate(p.begin(), p.begin() + 1, p.end());
        for(auto& p : pointsWithinGeometry)
            std::rotate(p.begin(), p.begin() + 1, p.end());
        for(auto& p : pointsOutsideGeometry)
            std::rotate(p.begin(), p.begin() + 1, p.end());

        checkAxisAlignedGeometry(cornerPoints, pointsWithinGeometry, pointsOutsideGeometry, i);
    }

    std::cout << "testing for non axis-aligned quadrilateral" << std::endl;

    cornerPoints[0][0] += 0.5;
    cornerPoints[1][0] -= 0.5;
    cornerPoints[2][0] += 0.5;
    cornerPoints[3][0] -= 0.5;

    GlobalPosition pointToCheck5 = {0.0, 0.5, 0.5};

    pointsWithinGeometry = {pointToCheck5};

    pointsOutsideGeometry[1] = pointsOutsideGeometry[0];

    permutatePointsAndTest(cornerPoints, pointsWithinGeometry, pointsOutsideGeometry);

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
