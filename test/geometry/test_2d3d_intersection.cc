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

    using Geometry3D = Dune::MultiLinearGeometry<double, 3, dimworld>;
    using Geometry2D = Dune::MultiLinearGeometry<double, 2, dimworld>;
    using Test = Dumux::GeometryIntersection<Geometry3D, Geometry2D>;
    typename Test::Intersection intersectionPolygon;

    std::vector<Dune::FieldVector<double, dimworld>> cubeCorners({
        {0.8619712261141845, 1.206565697102448, 0.4152254463380627},
        {0.8715547951379876, 1.175481516085758, 0.2599429674402531},
        {0.7046736652858085, 1.209207158782034, 0.4447521639960824},
        {0.7307137896935333, 1.171606500732276, 0.2574814772189353},
        {0.8926414294658755, 1,                 0.4012781039132332},
        {0.8316315415883838, 1,                 0.2581023823497354},
        {0.7349272235105835, 1,                 0.4613391568883559},
        {0.7479805294570721, 1,                 0.3155931875126768} });

    std::vector<Dune::FieldVector<double, dimworld>> quadCorners({
        {1.000000000000000000000000000000, 1.000000000000000000000000000000, 1.000000000000000000000000000000},
        {1.000000000000000000000000000000, 1.000000000000000000000000000000, 0.327819111366782434124900191819},
        {0.303996171086607036571081152942, 1.000000000000000000000000000000, 1.000000000000000000000000000000},
        {0.297170743844278550938042826601, 1.000000000000000000000000000000, 0.293625720570556747457402479995} });

    std::vector<Dune::FieldVector<double, dimworld>> quad2Corners({
        {1.000000000000000000000000000000, 1.000000000000000000000000000000, 0.327819111366782434124900191819},
        {1.000000000000000000000000000000, 1.000000000000000000000000000000, 0.0},
        {0.297170743844278550938042826601, 1.000000000000000000000000000000, 0.293625720570556747457402479995},
        {0.325413399309280815252520824288, 1.000000000000000000000000000000, 0.000000000000000000000000000000} });

    Geometry3D cube(Dune::GeometryTypes::cube(dimworld), cubeCorners);
    Geometry2D quad(Dune::GeometryTypes::cube(dimworld-1), quadCorners);

    // write out cube
    std::vector<std::vector<unsigned int>> faceCornerIndices({
        {0, 1, 3, 2},
        {0, 1, 5, 4},
        {1, 3, 7, 5},
        {0, 2, 6, 4},
        {2, 3, 7, 6},
        {4, 5, 7, 6}
    });
    for (unsigned int i = 0; i < 6; ++i)
    {
        std::vector<Point> corners;
        for (unsigned int c = 0; c < 4; ++c)
            corners.emplace_back(cubeCorners[faceCornerIndices[i][c]]);
        const auto t = Dumux::triangulate<2, dimworld>(corners);
        Dumux::writeVTKPolyDataTriangle(t, "cube_face_" + std::to_string(i));
    }

    // write out quadrilateral
    {
        std::vector<Point> quadRearrange({quadCorners[0], quadCorners[1], quadCorners[3], quadCorners[2]});
        const auto quadTriangulation = Dumux::triangulate<2, dimworld>(quadRearrange);
        Dumux::writeVTKPolyDataTriangle(quadTriangulation, "quad1");
    }
    {
        std::vector<Point> quadRearrange({quad2Corners[0], quad2Corners[1], quad2Corners[3], quad2Corners[2]});
        const auto quadTriangulation = Dumux::triangulate<2, dimworld>(quadRearrange);
        Dumux::writeVTKPolyDataTriangle(quadTriangulation, "quad2");
    }

    int numInts1 = 0;
    int numInts2 = 0;
    double area = 0.0;
    if (Test::intersection(cube, quad, intersectionPolygon))
    {
        const auto triangulation = Dumux::triangulate<2, dimworld>(intersectionPolygon);
        numInts1 = triangulation.size();

        // try to make a dune geometry with the result
        for (const auto& t : triangulation)
        {
            // might run into an assertion
            std::vector<Point> corners(t.begin(), t.end());
            Geometry2D triangle(Dune::GeometryTypes::simplex(dimworld-1), corners);
            std::cout << "Volume: " << triangle.volume() << std::endl;
            area += triangle.volume();
        }

        Dumux::writeVTKPolyDataTriangle(triangulation, "quad1_intersections");
    }

    quad = Geometry2D(Dune::GeometryTypes::cube(dimworld-1), quad2Corners);
    if (Test::intersection(cube, quad, intersectionPolygon))
    {
        const auto triangulation = Dumux::triangulate<2, dimworld>(intersectionPolygon);
        numInts2 = triangulation.size();

        // try to make a dune geometry with the result
        for (const auto& t : triangulation)
        {
            // might run into an assertion
            std::vector<Point> corners(t.begin(), t.end());
            Geometry2D triangle(Dune::GeometryTypes::simplex(dimworld-1), corners);
            std::cout << "Volume: " << triangle.volume() << std::endl;
            area += triangle.volume();
        }

        Dumux::writeVTKPolyDataTriangle(triangulation, "quad2_intersections");
    }

    // compute intersection face area
    std::vector<Dune::FieldVector<double, dimworld>> isFaceCorners({
    {0.8926414294658755, 1,                 0.4012781039132332},
    {0.8316315415883838, 1,                 0.2581023823497354},
    {0.7349272235105835, 1,                 0.4613391568883559},
    {0.7479805294570721, 1,                 0.3155931875126768} });
    Geometry2D isFace(Dune::GeometryTypes::cube(dimworld-1), isFaceCorners);
    const auto isFaceArea = isFace.volume();

    using std::abs;
    std::cout << "Intersection area mismatch: " << std::setprecision(30) << abs(area-isFaceArea) << std::endl;
    std::cout << "Num triangles 1: " << numInts1 << std::endl;
    std::cout << "Num triangles 2: " << numInts2 << std::endl;

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
