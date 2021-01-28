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
        {1.7837834300248889274, 1,                      0.49999996467938950806},
        {1.7835132521830019403, 1,                      0.54545451802928024421},
        {1.7662934121967939216, 0.9617398865062183333,  0.48484841728616778767},
        {1.7574133143618040354, 0.94260982975932738892, 0.49999992026450240656},
        {1.7935699809263809801, 0.97237560778921372506, 0.47177679079232959225},
        {1.7983281674561839569, 0.95856341168382064311, 0.48039248052374511344},
        {1.7780058298299370456, 0.95058662072157407152, 0.46746892371917841968},
        {1.7759899038176580888, 0.93411549429543205836, 0.47177676118240491343} });

    std::vector<Dune::FieldVector<double, dimworld>> quadCorners({
        {1.8749999999999489297, 1, 0.49999999999913791182},
        {1.9166666666666329899, 1, 0.58333333333253922781},
        {1.7499999999998980815, 1, 0.49999999999958189001},
        {1.8749999999999489297, 1, 0.62499999999946176388} });

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
    std::vector<Point> quadRearrange({quadCorners[0], quadCorners[1], quadCorners[3], quadCorners[2]});
    const auto quadTriangulation = Dumux::triangulate<2, dimworld>(quadRearrange);
    Dumux::writeVTKPolyDataTriangle(quadTriangulation, "quad");

    if (Test::intersection(cube, quad, intersectionPolygon))
    {
        const auto triangulation = Dumux::triangulate<2, dimworld>(intersectionPolygon);
        Dumux::writeVTKPolyDataTriangle(triangulation, "quad_intersections");
        DUNE_THROW(Dune::InvalidStateException, "Found unexpected intersection!");
    }

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
