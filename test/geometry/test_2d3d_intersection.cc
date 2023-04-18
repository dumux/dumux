//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// graham convex hull test + triangulation
#include <config.h>

#include <fstream>
#include <iostream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/geometry/geometryintersection.hh>
#include <dumux/geometry/triangulation.hh>
#include "writetriangulation.hh"

template<int dimworld = 3>
void testSegTriangle(const Dune::FieldVector<double, dimworld>& a,
                     const Dune::FieldVector<double, dimworld>& b,
                     const Dune::FieldVector<double, dimworld>& c,
                     const Dune::FieldVector<double, dimworld>& q,
                     const Dune::FieldVector<double, dimworld>& p,
                     bool expectIntersection = true)
{
    using GlobalPosition = Dune::FieldVector<double, dimworld>;
    using CornerStorage = std::vector<GlobalPosition>;
    using TriGeometry = Dune::MultiLinearGeometry<double, 2, dimworld>;
    using SegGeometry = Dune::MultiLinearGeometry<double, 1, dimworld>;
    using Policy = Dumux::IntersectionPolicy::PointPolicy<double, dimworld>;
    using Test = Dumux::GeometryIntersection<TriGeometry, SegGeometry, Policy>;
    typename Test::Intersection intersection;

    const auto tria = TriGeometry(Dune::GeometryTypes::triangle, CornerStorage({a, b, c}));
    const auto seg = SegGeometry(Dune::GeometryTypes::line, CornerStorage({q, p}));
    if (Test::intersection(tria, seg, intersection))
    {
        if (expectIntersection)
            std::cout << "Found intersection point: " << intersection << std::endl;
        else
            DUNE_THROW(Dune::InvalidStateException, "Found false positive: " << intersection);
    }
    else
    {
        if (expectIntersection)
            DUNE_THROW(Dune::InvalidStateException, "No intersection found!");
        else
            std::cout << "No intersection point (as expected). " << std::endl;
    }
}

template<typename ctype>
auto makeUnitCube(ctype scale = 1.0)
{
    using CubeGeometry = Dune::MultiLinearGeometry<ctype, 3, 3>;
    using GlobalPosition = typename CubeGeometry::GlobalCoordinate;

    const auto unit = scale*1.0;
    return CubeGeometry(
        Dune::GeometryTypes::cube(3),
        std::vector<GlobalPosition>({
            {0.0, 0.0, 0.0},  {unit, 0.0, 0.0},  {0.0, unit, 0.0},  {unit, unit, 0.0},
            {0.0, 0.0, unit}, {unit, 0.0, unit}, {0.0, unit, unit}, {unit, unit, unit}
        })
    );
}

template<typename ctype>
auto makeHorizontalQuad(ctype scale = 1.0, ctype zPositionUnscaled = 1.0)
{
    using QuadGeometry = Dune::MultiLinearGeometry<ctype, 2, 3>;
    using GlobalPosition = typename QuadGeometry::GlobalCoordinate;

    zPositionUnscaled *= scale;
    const auto unit = scale*1.0;
    return QuadGeometry(
        Dune::GeometryTypes::cube(2),
        std::vector<GlobalPosition>({
            {0.0, 0.0, zPositionUnscaled},  {unit, 0.0, zPositionUnscaled},
            {0.0, unit, zPositionUnscaled}, {unit, unit, zPositionUnscaled}
        })
    );
}

template<typename ctype>
auto makeInclinedTriangle(ctype scale = 1.0)
{
    using TriaGeometry = Dune::MultiLinearGeometry<ctype, 2, 3>;
    using GlobalPosition = typename TriaGeometry::GlobalCoordinate;
    return TriaGeometry(
        Dune::GeometryTypes::simplex(2),
        std::vector<GlobalPosition>({
            {-0.1*scale, -0.1*scale, 0.3*scale},
            {1.1*scale,  -0.1*scale, 0.3*scale},
            {0.5*scale,   2.0*scale, 0.8*scale}
        })
    );
}

int main(int argc, char* argv[])
{
    constexpr int dimworld = 3;
    using Point = Dune::FieldVector<double, dimworld>;

    for (auto scale : {1.0, 1e3, 1e12, 1e-12})
    {
        std::cout << "Testing scale " << scale << std::endl;

        const Point a{0.0*scale, 0.0*scale, 0.1*scale};
        const Point b{1.0*scale, 0.0*scale, 0.0*scale};
        const Point c{0.0*scale, 1.0*scale, 0.0*scale};

        const Point p0{0.5*scale, 0.5*scale, -1.0*scale};
        const Point q0{0.5*scale, 0.5*scale, 1.0*scale};

        testSegTriangle(a, b, c, p0, q0);
        testSegTriangle(a, b, c, q0, p0);

        const Point p1 = a;
        const Point q1 = b;

        testSegTriangle(a, b, c, p1, q1, false);
        testSegTriangle(a, b, c, q1, p1, false);

        const Point p2{0.0*scale, 0.0*scale, 0.0*scale};
        const Point q2{0.0*scale, 0.0*scale, 0.2*scale};

        testSegTriangle(a, b, c, p2, q2);

        const auto cube = makeUnitCube(scale);
        const auto quad = makeHorizontalQuad(scale, 0.5);
        using Geometry3D = decltype(cube);
        using Geometry2D = decltype(quad);
        using Test = Dumux::GeometryIntersection<Geometry3D, Geometry2D>;
        typename Test::Intersection intersectionPolygon;

        if (Test::intersection(cube, quad, intersectionPolygon))
        {
            const auto triangulation = Dumux::triangulate<2, dimworld>(intersectionPolygon);
            Dumux::writeVTKPolyDataTriangle(triangulation, "quad_intersections");
            if (triangulation.size() != 4)
                DUNE_THROW(Dune::InvalidStateException, "Found " << triangulation.size() << " instead of 4 intersections!");
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "No intersections found!");

        // quad that is slightly above the cube
        const auto close_quad = makeHorizontalQuad(scale, 1.0+1e-6);
        if (Test::intersection(cube, close_quad, intersectionPolygon))
            DUNE_THROW(Dune::InvalidStateException, "Found unexpected intersection with close quad");

        const auto tri = makeInclinedTriangle(scale);
        if (Test::intersection(cube, tri, intersectionPolygon))
        {
            const auto triangulation = Dumux::triangulate<2, dimworld>(intersectionPolygon);
            Dumux::writeVTKPolyDataTriangle(triangulation, "tri_intersections");
            if (triangulation.size() != 6)
                DUNE_THROW(Dune::InvalidStateException, "Found " << triangulation.size() << " instead of 6 intersections!");
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "No intersections found!");

        // Test a tricky situation intersecting a hexagon with two quadrilaterals.
        // One face of the hexagon is in-plane with the two quadrilaterals, which
        // together enclose the entire hexagon face. Thus, the sum of the areas of
        // the two intersections must be equal to that of the respective hexagon face.
        // The separating edge between the two quadrilaterals passes very close by one
        // of the hex face corners. This situation was taken from an actual mortar
        // application which exposed a bug in the 1d-2d algorithm in 3d space.
        std::vector<Dune::FieldVector<double, dimworld>> hexCorners({
            {0.8619712261141845*scale, 1.206565697102448*scale, 0.4152254463380627*scale},
            {0.8715547951379876*scale, 1.175481516085758*scale, 0.2599429674402531*scale},
            {0.7046736652858085*scale, 1.209207158782034*scale, 0.4447521639960824*scale},
            {0.7307137896935333*scale, 1.171606500732276*scale, 0.2574814772189353*scale},
            {0.8926414294658755*scale, 1.0*scale,               0.4012781039132332*scale},
            {0.8316315415883838*scale, 1.0*scale,               0.2581023823497354*scale},
            {0.7349272235105835*scale, 1.0*scale,               0.4613391568883559*scale},
            {0.7479805294570721*scale, 1.0*scale,               0.3155931875126768*scale}
        });

        std::vector<Dune::FieldVector<double, dimworld>> quad1Corners({
            {1.000000000000000000000000000000*scale, 1.0*scale, 1.000000000000000000000000000000*scale},
            {1.000000000000000000000000000000*scale, 1.0*scale, 0.327819111366782434124900191819*scale},
            {0.303996171086607036571081152942*scale, 1.0*scale, 1.000000000000000000000000000000*scale},
            {0.297170743844278550938042826601*scale, 1.0*scale, 0.293625720570556747457402479995*scale}
        });

        std::vector<Dune::FieldVector<double, dimworld>> quad2Corners({
            {1.000000000000000000000000000000*scale, 1.0*scale, 0.327819111366782434124900191819*scale},
            {1.000000000000000000000000000000*scale, 1.0*scale, 0.000000000000000000000000000000*scale},
            {0.297170743844278550938042826601*scale, 1.0*scale, 0.293625720570556747457402479995*scale},
            {0.325413399309280815252520824288*scale, 1.0*scale, 0.000000000000000000000000000000*scale}
        });

        Geometry3D hex(Dune::GeometryTypes::cube(dimworld), hexCorners);
        Geometry2D quad1(Dune::GeometryTypes::cube(dimworld-1), quad1Corners);
        Geometry2D quad2(Dune::GeometryTypes::cube(dimworld-1), quad2Corners);

        // lambda to compute the area of a triangulated intersection
        auto computeArea = [] (const auto& triangulation)
        {
            double a = 0.0;
            for (const auto& t : triangulation)
                a += Geometry2D(Dune::GeometryTypes::simplex(dimworld-1),
                                std::vector<Point>(t.begin(), t.end())).volume();
            return a;
        };

        double area = 0.0;
        if (Test::intersection(hex, quad1, intersectionPolygon))
        {
            const auto triangulation = Dumux::triangulate<2, dimworld>(intersectionPolygon);
            Dumux::writeVTKPolyDataTriangle(triangulation, "quad1_intersections");
            if (triangulation.size() != 5)
                DUNE_THROW(Dune::InvalidStateException, "Found " << triangulation.size() << " instead of 5 intersections!");

            area += computeArea(triangulation);
        }

        if (Test::intersection(hex, quad2, intersectionPolygon))
        {
            const auto triangulation = Dumux::triangulate<2, dimworld>(intersectionPolygon);
            Dumux::writeVTKPolyDataTriangle(triangulation, "quad2_intersections");
            if (triangulation.size() != 1)
                DUNE_THROW(Dune::InvalidStateException, "Found " << triangulation.size() << " instead of 1 intersections!");

            area += computeArea(triangulation);
        }

        // compute area of the intersecting hexagon face
        std::vector<Dune::FieldVector<double, dimworld>> hexFaceCorners({
            {0.8926414294658755*scale, 1.0*scale, 0.4012781039132332*scale},
            {0.8316315415883838*scale, 1.0*scale, 0.2581023823497354*scale},
            {0.7349272235105835*scale, 1.0*scale, 0.4613391568883559*scale},
            {0.7479805294570721*scale, 1.0*scale, 0.3155931875126768*scale}
        });
        Geometry2D hexFace(Dune::GeometryTypes::cube(dimworld-1), hexFaceCorners);
        const auto hexFaceArea = hexFace.volume();

        using std::abs;
        const auto areaMismatch = abs(area-hexFaceArea);
        std::cout << "Area mismatch is: " << areaMismatch << std::endl;
        if (areaMismatch > hexFaceArea*1e-7)
            DUNE_THROW(Dune::InvalidStateException, "Intersection area mismatch too large!");
    }

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
