#include <config.h>

#include <iostream>
#include <algorithm>
#include <initializer_list>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/geometry/geometryintersection.hh>

#ifndef DOXYGEN
Dune::MultiLinearGeometry<double, 1, 2>
makeLine(std::initializer_list<Dune::FieldVector<double, 2>>&& c, std::integral_constant<int, 2>)
{
    return {Dune::GeometryTypes::line, c};
}

Dune::MultiLinearGeometry<double, 1, 3>
makeLine(std::initializer_list<Dune::FieldVector<double, 2>>&& c, std::integral_constant<int, 3>)
{
    std::vector<Dune::FieldVector<double, 3>> corners;
    for (auto& corner : c)
    {
        Dune::FieldVector<double, 3> coord;
        coord[0] = corner[0];
        coord[1] = corner[1];
        coord[2] = 1.0;
        corners.emplace_back(std::move(coord));
    }

    return {Dune::GeometryTypes::line, corners};
}

template<int dimworld, class Policy>
bool testIntersection(const Dune::MultiLinearGeometry<double, 2, dimworld>& polygon,
                      const Dune::MultiLinearGeometry<double, 1, dimworld>& line,
                      bool foundExpected = true)
{
    using Test = Dumux::GeometryIntersection<Dune::MultiLinearGeometry<double, 2, dimworld>,
                                             Dune::MultiLinearGeometry<double, 1, dimworld>,
                                             Policy>;
    typename Test::Intersection intersection;
    bool found = Test::intersection(polygon, line, intersection);

    const std::string policyAppendix = Policy::dimIntersection == 0 ? "using point policy"
                                                                    : "using segment policy";

    if (!found && foundExpected)
        std::cerr << "Failed detecting intersection with " << line.corner(0) << " " << line.corner(1) << ", " << policyAppendix << std::endl;
    else if (found && foundExpected)
        std::cout << "Found intersection with " << line.corner(0) << " " << line.corner(1) << ", " << policyAppendix << std::endl;
    else if (found && !foundExpected)
        std::cerr << "Found false positive: intersection with " << line.corner(0) << " " << line.corner(1) << ", " << policyAppendix << std::endl;
    else if (!found && !foundExpected)
        std::cout << "No intersection with " << line.corner(0) << " " << line.corner(1) << ", " << policyAppendix << std::endl;
    return (found == foundExpected);
}

template<int dimworld, class Quadrilateral, class Triangle>
void performTests(std::vector<bool>& returns, const Quadrilateral& quad, const Triangle& triangle)
{
    // test with both point and segment policy
    using PointPolicy = Dumux::IntersectionPolicy::PointPolicy<double, dimworld>;
    using SegmentPolicy = Dumux::IntersectionPolicy::SegmentPolicy<double, dimworld>;

    std::integral_constant<int, dimworld> dm;

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {0.0, 0.0}, {1.0, 0.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {0.0, 0.0}, {1.0, 0.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {0.0, 0.0}, {0.0, 1.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {0.0, 0.0}, {0.0, 1.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {0.0, 0.0}, {1.0, 1.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {0.0, 0.0}, {1.0, 1.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {1.0, 0.0}, {1.0, 1.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {1.0, 0.0}, {1.0, 1.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {1.0, 1.0}, {0.0, 1.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {1.0, 1.0}, {0.0, 1.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {0.5, 0.0}, {0.5, 1.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {0.5, 0.0}, {0.5, 1.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {0.0, 0.5}, {1.0, 0.5} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {0.0, 0.5}, {1.0, 0.5} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {0.5, 0.5}, {0.5, 2.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {0.5, 0.5}, {0.5, 2.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {0.5, 0.5}, {0.5, -2.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {0.5, 0.5}, {0.5, -2.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {0.5, 0.5}, {-2.0, 0.5} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {0.5, 0.5}, {-2.0, 0.5} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {0.5, 0.5}, {2.0, 0.5} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {0.5, 0.5}, {2.0, 0.5} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {0.5, 0.0}, {0.5, -2.0} }, dm), false));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {0.5, 0.0}, {0.5, -2.0} }, dm)));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine({ {0.5, 1.0}, {0.5, 2.0} }, dm), false));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine({ {0.5, 1.0}, {0.5, 2.0} }, dm)));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine({ {0.0, 0.0}, {1.0, 0.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine({ {0.0, 0.0}, {1.0, 0.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine({ {0.0, 0.0}, {0.0, 1.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine({ {0.0, 0.0}, {0.0, 1.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine({ {0.0, 0.0}, {1.0, 1.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine({ {0.0, 0.0}, {1.0, 1.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine({ {0.5, 0.0}, {0.5, 1.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine({ {0.5, 0.0}, {0.5, 1.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine({ {0.0, 0.5}, {1.0, 0.5} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine({ {0.0, 0.5}, {1.0, 0.5} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine({ {0.5, 0.5}, {0.0, 0.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine({ {0.5, 0.5}, {0.0, 0.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine({ {0.0, 0.0}, {0.5, 0.5} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine({ {0.0, 0.0}, {0.5, 0.5} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine({ {0.5, 0.5}, {-2.0, 0.5} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine({ {0.5, 0.5}, {-2.0, 0.5} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine({ {0.5, 0.5}, {0.5, -2.0} }, dm)));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine({ {0.5, 0.5}, {0.5, -2.0} }, dm), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine({ {1.0, 1.0}, {0.0, 1.0} }, dm), false));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine({ {1.0, 1.0}, {0.0, 1.0} }, dm)));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine({ {1.0, 0.0}, {1.0, 1.0} }, dm), false));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine({ {1.0, 0.0}, {1.0, 1.0} }, dm)));
}

#endif

int main (int argc, char *argv[])
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // collect returns to determine exit code
    std::vector<bool> returns;

    // tests with dimWorld = 2
    {
        constexpr int dimworld = 2;
        constexpr int dim = 2;

        // we test quadrilateral-line & triangle-line intersections
        using CornerStorage = std::vector<Dune::FieldVector<double, dimworld>>;
        CornerStorage quadCorners({ {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0} });
        CornerStorage triaCorners({ {0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0} });

        using Geometry = Dune::MultiLinearGeometry<double, dim, dimworld>;
        Geometry quad(Dune::GeometryTypes::cube(dim), quadCorners);
        Geometry triangle(Dune::GeometryTypes::simplex(dim), triaCorners);

        std::cout << "Testing algorithms for dimworld = 2" << std::endl;
        performTests<dimworld>(returns, quad, triangle);
    }

    // tests with dimWorld = 3
    {
        constexpr int dimworld = 3;
        constexpr int dim = 2;

        // we test quadrilateral-line & triangle-line intersections
        using CornerStorage = std::vector<Dune::FieldVector<double, dimworld>>;
        CornerStorage quadCorners({ {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0} });
        CornerStorage triaCorners({ {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0} });

        using Geometry = Dune::MultiLinearGeometry<double, dim, dimworld>;
        Geometry quad(Dune::GeometryTypes::cube(dim), quadCorners);
        Geometry triangle(Dune::GeometryTypes::simplex(dim), triaCorners);

        std::cout << "Testing algorithms for dimworld = 3" << std::endl;
        performTests<dimworld>(returns, quad, triangle);
    }

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), [](bool i){ return !i; }))
        return 1;

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
