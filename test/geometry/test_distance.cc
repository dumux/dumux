// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for distance computations.
 */
#include <config.h>

#include <array>
#include <vector>
#include <random>
#include <cmath>

#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/geometry/circumsphere.hh>
#include <dumux/geometry/distance.hh>
#include <dumux/geometry/distancefield.hh>
#include <dumux/geometry/normal.hh>
#include <dumux/io/format.hh>

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

// helper function to make a triangle geometry from field vectors
template<class Triangle, int dimWorld>
Triangle makeTriangle(const Dune::FieldVector<typename Triangle::ctype, dimWorld>& p1,
                      const Dune::FieldVector<typename Triangle::ctype, dimWorld>& p2,
                      const Dune::FieldVector<typename Triangle::ctype, dimWorld>& p3)
{
    using GlobalPosition = Dune::FieldVector<typename Triangle::ctype, dimWorld>;
    using Corners = std::vector<GlobalPosition>;
    return { Dune::GeometryTypes::simplex(2), Corners{{p1, p2, p3}} };
}

// helper function to make a quadrilateral geometry from field vectors
template<class Quadrilateral, int dimWorld>
Quadrilateral makeQuadrilateral(const Dune::FieldVector<typename Quadrilateral::ctype, dimWorld>& p1,
                                const Dune::FieldVector<typename Quadrilateral::ctype, dimWorld>& p2,
                                const Dune::FieldVector<typename Quadrilateral::ctype, dimWorld>& p3,
                                const Dune::FieldVector<typename Quadrilateral::ctype, dimWorld>& p4)
{
    using GlobalPosition = Dune::FieldVector<typename Quadrilateral::ctype, dimWorld>;
    using Corners = std::vector<GlobalPosition>;
    return { Dune::GeometryTypes::quadrilateral, Corners{{p1, p2, p3, p4}} };
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

template<class Geometry>
void checkDistanceField(const std::vector<Geometry>& geometries,
                        const typename Geometry::GlobalCoordinate& p,
                        typename Geometry::ctype expectedDistance,
                        int expectedIndex)
{
    Dumux::AABBDistanceField<Geometry> distanceField(geometries);
    const auto [d, idx] = distanceField.distanceAndIndex(p);

    if (Dune::FloatCmp::ne(d, expectedDistance))
        DUNE_THROW(Dune::InvalidStateException,
            Dumux::Fmt::format(
                "Distance field returned wrong minimal distance. Expected {}, got {}",
                expectedDistance, d
            )
        );

    if (idx != expectedIndex)
        DUNE_THROW(Dune::InvalidStateException,
            Dumux::Fmt::format(
                "Distance field returned wrong index for minimal distance. Expected {}, got {}",
                expectedIndex, idx
            )
        );
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

            // check distance field class (2D), last segment is closest to p
            const auto segment4 = makeSegment<Segment>(origin.corner(0) + 0.5*n, p.corner(0) + 0.5*n);
            const std::vector<Segment> segments{segment, segment2, segment3, segment4};
            checkDistanceField(segments, origin.corner(0), 0.5*scale, segments.size() - 1);

            // check distance point-triangle
            if constexpr(GlobalPosition::dimension == 3)
            {
                // create a triangle whose normal is aligned with (p-origin), touching the sphere
                using Points = std::array<GlobalPosition, 3>;
                auto rotationAxis = p.corner(0); rotationAxis /= rotationAxis.two_norm();
                const auto rotate1 = make3DTransformation(1.0, origin.corner(0), rotationAxis, 2.0/3.0*M_PI, /*verbose*/false);
                const auto rotate2 = make3DTransformation(1.0, origin.corner(0), rotationAxis, 4.0/3.0*M_PI, /*verbose*/false);
                const auto firstTrianglePoints = Points{
                    p.corner(0) + n, rotate1(n) + p.corner(0), rotate2(n) + p.corner(0)
                };

                checkGeometryDistance(
                    scale,
                    distancePointTriangle(
                        origin.corner(0), firstTrianglePoints[0], firstTrianglePoints[1], firstTrianglePoints[2]
                    ),
                    "point-triangle"
                );

                // check convenience function
                using Polygon = OtherGeometryType<double, 2, coorddim>;
                const auto firstTriangle = makeTriangle<Polygon, coorddim>(
                    firstTrianglePoints[0], firstTrianglePoints[1], firstTrianglePoints[2]
                );

                checkGeometryDistance(scale, distance(origin, firstTriangle), "point-triangle");
                checkGeometryDistance(scale, distance(firstTriangle, origin), "triangle-point");

                // check circumsphere function
                const auto sphere = circumSphereOfTriangle(firstTriangle);
                const auto& expectedCenter = firstTriangle.center();
                const auto expectedRadius = (firstTriangle.corner(0) - expectedCenter).two_norm();
                if ((sphere.center() - expectedCenter).two_norm() > 1e-12*expectedRadius)
                    DUNE_THROW(Dune::InvalidStateException,
                        "Wrong circumsphere center. Expected "
                        << expectedCenter << " got " << sphere.center()
                    );

                if (Dune::FloatCmp::ne(sphere.radius(), expectedRadius))
                    DUNE_THROW(Dune::InvalidStateException, Fmt::format(
                        "Wrong circumsphere radius Expected {}, got {}",
                        expectedRadius, sphere.radius()
                    ));

                // check quadrilateral
                const auto& [a, b, c] = firstTrianglePoints;
                const auto d = b + (b - a);
                const auto quadrilateral = makeQuadrilateral<Polygon, coorddim>(a, b, c, d);
                checkGeometryDistance(scale, distance(origin, quadrilateral), "point-quadrilateral");
                checkGeometryDistance(scale, distance(quadrilateral, origin), "quadrilateral-point");

                // shift the first triangle
                const auto secondTrianglePoints = Points{
                    segment.corner(0), segment.corner(1), 0.5*(segment.corner(0) + segment.corner(1)) + n
                };
                const auto secondTriangle = makeTriangle<Polygon, coorddim>(
                    secondTrianglePoints[0], secondTrianglePoints[1], secondTrianglePoints[2]
                );
                checkGeometryDistance(scale, distancePointTriangle(origin.corner(0), secondTriangle), "point-triangle");

                // create a triangle not touching the sphere
                const auto thirdTriangle = Points{
                    segment2.corner(0), segment2.corner(1), 0.5*(segment2.corner(0) + segment2.corner(1)) + n
                };

                checkGeometryDistance(
                    sqrt(2.0)*scale,
                    distancePointTriangle(
                        p.corner(0), thirdTriangle[0], thirdTriangle[1], thirdTriangle[2]
                    ),
                    "point-triangle"
                );

                // check the distance field class
                const auto closestQuadrilateral = makeQuadrilateral<Polygon, coorddim>(0.5*a, 0.5*b, 0.5*c, 0.5*d);
                std::vector<Polygon> polygons{firstTriangle, quadrilateral, secondTriangle, closestQuadrilateral};
                checkDistanceField(polygons, origin.corner(0), 0.5*scale, polygons.size() - 1);
            }
        }
    }
}

// checks the average distance of a line and of segments on the z-axis to a box
template<class Hexahedron>
void runAverageDistanceTests()
{
    using namespace Dumux;
    using GlobalPosition = typename Hexahedron::GlobalCoordinate;

    // the box [1, 2] x [0, 1] x [0, 2] next to the z-axis
    std::vector<GlobalPosition> corners;
    for (double z : {0.0, 2.0})
        for (double y : {0.0, 1.0})
            for (double x : {1.0, 2.0})
                corners.push_back({x, y, z});
    const Hexahedron box(Dune::GeometryTypes::cube(3), corners);

    // the distance to the z-axis is r = sqrt(x^2 + y^2) for all z, so the exact average over the box
    // is the average of r over the square [1, 2] x [0, 1], computed with the antiderivative F of r
    const auto F = [](double x, double y){
        using std::hypot; using std::log;
        const auto r = hypot(x, y);
        return (2.0*x*y*r + x*x*x*log(y + r) + y*y*y*log(x + r))/6.0;
    };
    const double expected = F(2.0, 1.0) - F(1.0, 1.0) - F(2.0, 0.0) + F(1.0, 0.0);

    // segments on the z-axis: covering the box height, covering a part of it, and above the box
    const auto onAxis = [](double z){ return GlobalPosition({0.0, 0.0, z}); };
    const auto covering = std::array{onAxis(0.0), onAxis(2.0)};
    const auto partial = std::array{onAxis(0.5), onAxis(1.0)};
    const auto above = std::array{onAxis(3.0), onAxis(4.0)};

    // the line distance does not depend on the segment: it matches the exact value with a high
    // integration order (r is smooth on the box) and is identical for all segments with the default order
    const auto lineDistance = averageDistanceLineGeometry(covering[0], covering[1], box);
    for (const auto& [a, b] : {covering, partial, above})
    {
        const auto computed = averageDistanceLineGeometry(a, b, box, 15);
        if (Dune::FloatCmp::ne(expected, computed, 1e-12))
            DUNE_THROW(Dune::InvalidStateException, Fmt::format(
                "Unexpected average line-box distance. Expected {}, got {}", expected, computed
            ));
        checkGeometryDistance(lineDistance, averageDistanceLineGeometry(a, b, box), "average line-box");
    }

    // the segment distance equals the line distance only if the segment covers the box height,
    // otherwise the parts of the box beyond the segment ends are farther away
    checkGeometryDistance(lineDistance, averageDistanceSegmentGeometry(covering[0], covering[1], box), "average segment-box");
    for (const auto& [a, b] : {partial, above})
        if (!(averageDistanceSegmentGeometry(a, b, box) > lineDistance))
            DUNE_THROW(Dune::InvalidStateException,
                "Average segment-box distance has to exceed the line distance if the segment does not cover the box height");
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

    // average distances of a line and of segments to a box
    runAverageDistanceTests<Dune::MultiLinearGeometry<double, 3, 3>>();

    return 0;
}
