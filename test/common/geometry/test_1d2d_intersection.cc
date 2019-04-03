#include <config.h>

#include <iostream>
#include <algorithm>
#include <initializer_list>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/geometry/geometryintersection.hh>


#ifndef DOXYGEN
template<int dimworld>
Dune::MultiLinearGeometry<double, 1, dimworld>
makeLine(std::initializer_list<Dune::FieldVector<double, dimworld>>&& c,
         double thirdCoord = 1.0)
{
    if (dimworld == 2)
        return {Dune::GeometryTypes::line, c};

    // add value for the third coordinate
    else if (dimworld == 3)
    {
        std::vector<Dune::FieldVector<double, dimworld>> corners;
        for (auto& corner : c)
        {
            Dune::FieldVector<double, dimworld> coord;
            coord[0] = corner[0];
            coord[1] = corner[1];
            coord[2] = thirdCoord;
            corners.emplace_back(std::move(coord));
        }

        return {Dune::GeometryTypes::line, corners};
    }
    else
        DUNE_THROW(Dune::InvalidStateException, "dimWorld must be greater than 2 or 3");
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

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {0.0, 0.0}, {1.0, 0.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {0.0, 0.0}, {1.0, 0.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {0.0, 0.0}, {0.0, 1.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {0.0, 0.0}, {0.0, 1.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {0.0, 0.0}, {1.0, 1.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {0.0, 0.0}, {1.0, 1.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {1.0, 0.0}, {1.0, 1.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {1.0, 0.0}, {1.0, 1.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {1.0, 1.0}, {0.0, 1.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {1.0, 1.0}, {0.0, 1.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {0.5, 0.0}, {0.5, 1.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {0.5, 0.0}, {0.5, 1.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {0.0, 0.5}, {1.0, 0.5} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {0.0, 0.5}, {1.0, 0.5} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {0.5, 0.5}, {0.5, 2.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {0.5, 0.5}, {0.5, 2.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {0.5, 0.5}, {0.5, -2.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {0.5, 0.5}, {0.5, -2.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {0.5, 0.5}, {-2.0, 0.5} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {0.5, 0.5}, {-2.0, 0.5} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {0.5, 0.5}, {2.0, 0.5} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {0.5, 0.5}, {2.0, 0.5} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {0.5, 0.0}, {0.5, -2.0} }), false));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {0.5, 0.0}, {0.5, -2.0} })));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(quad, makeLine<dimworld>({ {0.5, 1.0}, {0.5, 2.0} }), false));
    returns.push_back(testIntersection<dimworld, PointPolicy>(quad, makeLine<dimworld>({ {0.5, 1.0}, {0.5, 2.0} })));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine<dimworld>({ {0.0, 0.0}, {1.0, 0.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine<dimworld>({ {0.0, 0.0}, {1.0, 0.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine<dimworld>({ {0.0, 0.0}, {0.0, 1.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine<dimworld>({ {0.0, 0.0}, {0.0, 1.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine<dimworld>({ {0.0, 0.0}, {1.0, 1.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine<dimworld>({ {0.0, 0.0}, {1.0, 1.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine<dimworld>({ {0.5, 0.0}, {0.5, 1.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine<dimworld>({ {0.5, 0.0}, {0.5, 1.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine<dimworld>({ {0.0, 0.5}, {1.0, 0.5} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine<dimworld>({ {0.0, 0.5}, {1.0, 0.5} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine<dimworld>({ {0.5, 0.5}, {0.0, 0.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine<dimworld>({ {0.5, 0.5}, {0.0, 0.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine<dimworld>({ {0.0, 0.0}, {0.5, 0.5} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine<dimworld>({ {0.0, 0.0}, {0.5, 0.5} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine<dimworld>({ {0.5, 0.5}, {-2.0, 0.5} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine<dimworld>({ {0.5, 0.5}, {-2.0, 0.5} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine<dimworld>({ {0.5, 0.5}, {0.5, -2.0} })));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine<dimworld>({ {0.5, 0.5}, {0.5, -2.0} }), false));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine<dimworld>({ {1.0, 1.0}, {0.0, 1.0} }), false));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine<dimworld>({ {1.0, 1.0}, {0.0, 1.0} })));

    returns.push_back(testIntersection<dimworld, SegmentPolicy>(triangle, makeLine<dimworld>({ {1.0, 0.0}, {1.0, 1.0} }), false));
    returns.push_back(testIntersection<dimworld, PointPolicy>(triangle, makeLine<dimworld>({ {1.0, 0.0}, {1.0, 1.0} })));
}

#endif

int main (int argc, char *argv[]) try
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
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
