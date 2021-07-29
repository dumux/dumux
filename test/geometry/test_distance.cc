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
#include "transformation.hh"

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
template<template<class, int, int> class PointType,
         template<class, int, int> class OtherGeometryType,
         int coorddim>
void runTests()
{
    using namespace Dumux;

    using Point = PointType<double, 0, coorddim>;
    using Segment = OtherGeometryType<double, 1, coorddim>;
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
            checkGeometryDistance(scale, distance(segment, p), "segment-point");

            // test point-segment distance (where projection is NOT on segment)
            auto v = p.corner(0);

            using std::sqrt;
            auto segment2 = makeSegment<Segment>(segment.corner(0) - v, segment.corner(1) - v);
            checkGeometryDistance(sqrt(2.0)*scale, distance(segment2, p), "segment-point");

            v *= 2.0;
            auto segment3 = makeSegment<Segment>(segment.corner(0) + v, segment.corner(1) + v);
            checkGeometryDistance(sqrt(2.0)*scale, distance(segment3, p), "segment-point");

            // for lines, we should always get a distance of scale
            checkGeometryDistance(scale, distancePointLine(p.corner(0), segment.corner(0), segment.corner(1)), "point-line");
            checkGeometryDistance(scale, distancePointLine(p.corner(0), segment2.corner(0), segment2.corner(1)), "point-line");
            checkGeometryDistance(scale, distancePointLine(p.corner(0), segment3.corner(0), segment3.corner(1)), "point-line");

            // check distance point-triangle
            if constexpr(GlobalPosition::dimension == 3)
            {
                // create a triangle whose normal is aligned with (p-origin), touching the sphere
                using Points = std::array<GlobalPosition, 3>;
                auto rotationAxis = p.corner(0); rotationAxis /= rotationAxis.two_norm();
                const auto rotate1 = make3DTransformation(1.0, origin.corner(0), rotationAxis, 2.0/3.0*M_PI, /*verbose*/false);
                const auto rotate2 = make3DTransformation(1.0, origin.corner(0), rotationAxis, 4.0/3.0*M_PI, /*verbose*/false);
                const auto firstTriangle = Points{p.corner(0) + n,
                                                  rotate1(n) + p.corner(0),
                                                  rotate2(n) + p.corner(0)};
                checkGeometryDistance(scale, distancePointTriangle(origin.corner(0), firstTriangle[0], firstTriangle[1], firstTriangle[2]), "point-triangle");

                // check convenience function
                using Polygon = OtherGeometryType<double, 2, coorddim>;
                const auto triangleGeometry = Polygon(Dune::GeometryTypes::simplex(2), std::vector<GlobalPosition>{firstTriangle.begin(), firstTriangle.end()});
                checkGeometryDistance(scale, distance(origin, triangleGeometry), "point-triangle");
                checkGeometryDistance(scale, distance(triangleGeometry, origin), "point-triangle");

                // check quadrilateral
                const auto& a = firstTriangle[0];
                const auto& b = firstTriangle[1];
                const auto& c = firstTriangle[2];
                const auto d = b + (b - a);
                const auto quadGeometry = Polygon(Dune::GeometryTypes::quadrilateral, std::vector<GlobalPosition>{a, b, c, d});
                checkGeometryDistance(scale, distance(origin, quadGeometry), "point-quadrilateral");
                checkGeometryDistance(scale, distance(quadGeometry, origin), "point-quadrilateral");

                // shift the first triangle
                const auto secondTriangle = Points{segment.corner(0),
                                                   segment.corner(1),
                                                   0.5*(segment.corner(0) + segment.corner(1)) + n};
                checkGeometryDistance(scale, distancePointTriangle(origin.corner(0), secondTriangle[0], secondTriangle[1], secondTriangle[2]), "point-triangle");
                checkGeometryDistance(scale, distancePointTriangle(p.corner(0), secondTriangle[0], secondTriangle[1], secondTriangle[2]), "point-triangle");

                // create a triangle not touching the sphere
                const auto thirdTriangle = Points{segment2.corner(0),
                                                  segment2.corner(1),
                                                  0.5*(segment2.corner(0) + segment2.corner(1)) + n};
                checkGeometryDistance(sqrt(2.0)*scale, distancePointTriangle(p.corner(0), thirdTriangle[0], thirdTriangle[1], thirdTriangle[2]), "point-triangle");
            }
        }
    }
}

// alias to ignore fourth optional template argument of Dune::MultiLinearGeometry
template<class ct, int mydim, int codim>
using DuneMultiLinearGeometry = Dune::MultiLinearGeometry<ct, mydim, codim>;

int main(int argc, char** argv)
{
    // affine geometries (2d)
    runTests<Dune::AffineGeometry, Dune::AffineGeometry, 2>();

    // affine geometries (3d)
    runTests<Dune::AffineGeometry, Dune::AffineGeometry, 3>();

    // multilinear geometries (2d)
    runTests<DuneMultiLinearGeometry, DuneMultiLinearGeometry, 2>();

    // multilinear geometries (3d)
    runTests<DuneMultiLinearGeometry, DuneMultiLinearGeometry, 3>();

    // mixed types (2d)
    runTests< Dune::AffineGeometry, DuneMultiLinearGeometry, 2>();
    runTests< DuneMultiLinearGeometry, Dune::AffineGeometry, 2>();

    // mixed types (3d)
    runTests< Dune::AffineGeometry, DuneMultiLinearGeometry, 3>();
    runTests< DuneMultiLinearGeometry, Dune::AffineGeometry, 3>();

    return 0;
}
