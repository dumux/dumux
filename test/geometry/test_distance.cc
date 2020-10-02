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
 * \brief Test for distance computations.
 */
#include <config.h>

#include <vector>
#include <random>
#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/geometry/distance.hh>
#include <dumux/geometry/normal.hh>

// helper function to make point geometry from field vector
template<class Point, int dimWorld>
Point makePoint(const Dune::FieldVector<typename Point::ctype, dimWorld>& p)
{
    using GlobalPosition = Dune::FieldVector<typename Point::ctype, dimWorld>;
    using Corners = std::vector<GlobalPosition>;
    return { Dune::GeometryTypes::vertex, Corners{{p}} };
}

// helper function to make a segment geometry from field vectors
template<class Segment, int dimWorld>
Segment makeSegment(const Dune::FieldVector<typename Segment::ctype, dimWorld>& p1,
                    const Dune::FieldVector<typename Segment::ctype, dimWorld>& p2)
{
    using GlobalPosition = Dune::FieldVector<typename Segment::ctype, dimWorld>;
    using Corners = std::vector<GlobalPosition>;
    return { Dune::GeometryTypes::line, Corners{{p1, p2}} };
}

// sample a point on a sphere with the given radius
template<class Point>
Point samplePointOnSphere(typename Point::ctype radius)
{
    static std::default_random_engine generator(std::random_device{}());
    static constexpr int dimWorld = Point::coorddimension;
    static_assert(dimWorld > 1, "This only works in 2d or 3d");

    // sample a point from random coordinates
    std::uniform_real_distribution<double> distro(-radius, radius);

    auto p = dimWorld == 2 ? makePoint<Point, dimWorld>({distro(generator), distro(generator)})
                           : makePoint<Point, dimWorld>({distro(generator), distro(generator), distro(generator)});

    // project it onto the sphere
    auto pos = p.corner(0);
    pos *= radius/pos.two_norm();
    p = makePoint<Point>(pos);

    return p;
}

// checks the result of a distance computation
template<class ctype>
void checkGeometryDistance(ctype expected, ctype computed, const std::string& geometryPairName)
{
    if (Dune::FloatCmp::ne(expected, computed))
        DUNE_THROW(Dune::InvalidStateException, "Unexpected " << geometryPairName << " distance ");
}

// checks the distances between various points with points/segments/lines
template<class Point, class Segment>
void runTests()
{
    using namespace Dumux;

    using ctype = typename Point::ctype;
    using GlobalPosition = typename Point::GlobalCoordinate;

    const auto origin = makePoint<Point>(GlobalPosition(0.0));
    for (ctype scale : {1e-12, 1.0, 1e12})
    {
        // check distances for some random points in a known distance
        for (int i = 0; i < 20; ++i)
        {
            const auto p = samplePointOnSphere<Point>(scale);

            // test point-point distance
            checkGeometryDistance(scale, distance(origin, p), "point-point");

            // test point-segment distance (where projection is on the segment)
            auto n = normal(p.corner(0));
            n *= scale/n.two_norm();

            auto segment = makeSegment<Segment>(origin.corner(0) + n, p.corner(0) + n);
            checkGeometryDistance(scale, distance(segment, p), "point-segment");

            // test point-segment distance (where projection is NOT on segment)
            auto v = p.corner(0);

            using std::sqrt;
            auto segment2 = makeSegment<Segment>(segment.corner(0) - v, segment.corner(1) - v);
            checkGeometryDistance(sqrt(2.0)*scale, distance(segment2, p), "segment-segment");

            v *= 2.0;
            auto segment3 = makeSegment<Segment>(segment.corner(0) + v, segment.corner(1) + v);
            checkGeometryDistance(sqrt(2.0)*scale, distance(segment2, p), "segment-segment");

            // for lines, we should always get a distance of scale
            checkGeometryDistance(scale, distancePointLine(p.corner(0), segment.corner(0), segment.corner(1)), "segment-segment");
            checkGeometryDistance(scale, distancePointLine(p.corner(0), segment2.corner(0), segment2.corner(1)), "segment-segment");
            checkGeometryDistance(scale, distancePointLine(p.corner(0), segment3.corner(0), segment3.corner(1)), "segment-segment");
        }
    }
}

int main(int argc, char** argv)
{
    // affine geometries (2d)
    runTests< Dune::AffineGeometry<double, 0, 2>, Dune::AffineGeometry<double, 1, 2> >();

    // affine geometries (3d)
    runTests< Dune::AffineGeometry<double, 0, 3>, Dune::AffineGeometry<double, 1, 3> >();

    // multilinear geometries (2d)
    runTests< Dune::MultiLinearGeometry<double, 0, 2>, Dune::MultiLinearGeometry<double, 1, 2> >();

    // multilinear geometries (3d)
    runTests< Dune::MultiLinearGeometry<double, 0, 3>, Dune::MultiLinearGeometry<double, 1, 3> >();

    // mixed types (2d)
    runTests< Dune::AffineGeometry<double, 0, 2>, Dune::MultiLinearGeometry<double, 1, 2> >();
    runTests< Dune::MultiLinearGeometry<double, 0, 2>, Dune::AffineGeometry<double, 1, 2> >();

    // mixed types (3d)
    runTests< Dune::AffineGeometry<double, 0, 3>, Dune::MultiLinearGeometry<double, 1, 3> >();
    runTests< Dune::MultiLinearGeometry<double, 0, 3>, Dune::AffineGeometry<double, 1, 3> >();

    return 0;
}
