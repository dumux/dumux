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

int main(int argc, char* argv[])
{
    constexpr int dimworld = 3;
    using Point = Dune::FieldVector<double, dimworld>;

    Point a{0.0, 0.0, 0.1};
    Point b{1.0, 0.0, 0.0};
    Point c{0.0, 1.0, 0.0};

    Point p0{0.5, 0.5, -1.0};
    Point q0{0.5, 0.5, 1.0};

    testSegTriangle(a, b, c, p0, q0);
    testSegTriangle(a, b, c, q0, p0);

    Point p1 = a;
    Point q1 = b;

    testSegTriangle(a, b, c, p1, q1, false);
    testSegTriangle(a, b, c, q1, p1, false);

    Point p2{0.0, 0.0, 0.0};
    Point q2{0.0, 0.0, 0.2};

    testSegTriangle(a, b, c, p2, q2);

    using Geometry3D = Dune::MultiLinearGeometry<double, 3, dimworld>;
    using Geometry2D = Dune::MultiLinearGeometry<double, 2, dimworld>;
    using Test = Dumux::GeometryIntersection<Geometry3D, Geometry2D>;
    typename Test::Intersection intersectionPolygon;

    std::vector<Dune::FieldVector<double, dimworld>> cubeCorners({
        {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0},
        {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}
    });

    std::vector<Dune::FieldVector<double, dimworld>> quadCorners({
        {0.0, 0.0, 0.5}, {1.0, 0.0, 0.5}, {0.0, 1.0, 0.5}, {1.0, 1.0, 0.5}
    });

    std::vector<Dune::FieldVector<double, dimworld>> triCorners({
        {-0.1, -0.1, 0.3}, {1.1, -0.1, 0.3}, {0.5, 2.0, 0.8}
    });

    Geometry3D cube(Dune::GeometryTypes::cube(dimworld), cubeCorners);
    Geometry2D quad(Dune::GeometryTypes::cube(dimworld-1), quadCorners);
    Geometry2D tri(Dune::GeometryTypes::simplex(dimworld-1), triCorners);

    if (Test::intersection(cube, quad, intersectionPolygon))
    {
        const auto triangulation = Dumux::triangulate<2, dimworld>(intersectionPolygon);
        Dumux::writeVTKPolyDataTriangle(triangulation, "quad_intersections");
        if (triangulation.size() != 4)
            DUNE_THROW(Dune::InvalidStateException, "Found " << triangulation.size() << " instead of 4 intersections!");
    }
    else
        DUNE_THROW(Dune::InvalidStateException, "No intersections found!");

    if (Test::intersection(cube, tri, intersectionPolygon))
    {
        const auto triangulation = Dumux::triangulate<2, dimworld>(intersectionPolygon);
        Dumux::writeVTKPolyDataTriangle(triangulation, "tri_intersections");
        if (triangulation.size() != 6)
            DUNE_THROW(Dune::InvalidStateException, "Found " << triangulation.size() << " instead of 6 intersections!");
    }
    else
        DUNE_THROW(Dune::InvalidStateException, "No intersections found!");

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
