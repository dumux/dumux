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
        {1.7837834300248889274342900535,  1,                                0.499999964679389508059870195211},
        {1.78351325218300194030973671033, 1,                                0.545454518029280244206802308327},
        {1.7662934121967939216091281196,  0.961739886506218333295237243874, 0.484848417286167787665362993721},
        {1.75741331436180403535729510622, 0.942609829759327388920553403295, 0.499999920264502406563877912049},
        {1.79356998092638098007967073499, 0.972375607789213725062893445283, 0.471776790792329592250098357908},
        {1.79832816745618395692929425422, 0.958563411683820643105491399183, 0.48039248052374511344098095833},
        {1.77800582982993704561636150174, 0.950586620721574071524173632497, 0.467468923719178419684538994261},
        {1.77598990381765808876934897853, 0.934115494295432058358130689157, 0.471776761182404913430588067058} });

    std::vector<Dune::FieldVector<double, dimworld>> quadCorners({
        {1.8749999999999489297408672428,  1, 0.499999999999137911821378565946},
        {1.91666666666663298990158637025, 1, 0.583333333332539227811253113032},
        {1.74999999999989808152633941063, 1, 0.499999999999581890008926166047},
        {1.8749999999999489297408672428,  1, 0.624999999999461763877661724109} });

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

        // try to make a dune geometry with the result
        for (const auto& t : triangulation)
        {
            // might run into an assertion
            std::vector<Point> corners(t.begin(), t.end());
            Geometry2D triangle(Dune::GeometryTypes::simplex(dimworld-1), corners);
            std::cout << "Volume: " << triangle.volume() << std::endl;
        }

        Dumux::writeVTKPolyDataTriangle(triangulation, "quad_intersections");
        DUNE_THROW(Dune::InvalidStateException, "Found unexpected intersection!");
    }

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
