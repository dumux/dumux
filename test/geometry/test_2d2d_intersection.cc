//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/geometry/geometryintersection.hh>

#ifndef DOXYGEN
template<int dimWorld>
Dune::FieldVector<double, dimWorld>
makePosition( double x, double y )
{
    if (dimWorld == 2)
        return {x, y};
    else if (dimWorld == 3)
        return {x, y, 1.0};
    else
        DUNE_THROW(Dune::InvalidStateException, "Invalid dimworld");
}

template<int dimWorld>
Dune::MultiLinearGeometry<double, 2, dimWorld>
makeTriangle( Dune::FieldVector<double, dimWorld>&& a,
              Dune::FieldVector<double, dimWorld>&& b,
              Dune::FieldVector<double, dimWorld>&& c )
{
    using CornerStorage = std::vector<Dune::FieldVector<double, dimWorld>>;
    return { Dune::GeometryTypes::triangle,
             CornerStorage({std::move(a), std::move(b), std::move(c)}) };
}

template<int dimWorld>
Dune::MultiLinearGeometry<double, 2, dimWorld>
makeQuadrilateral( Dune::FieldVector<double, dimWorld>&& a,
                   Dune::FieldVector<double, dimWorld>&& b,
                   Dune::FieldVector<double, dimWorld>&& c,
                   Dune::FieldVector<double, dimWorld>&& d )
{
    using CornerStorage = std::vector<Dune::FieldVector<double, dimWorld>>;
    return { Dune::GeometryTypes::quadrilateral,
             CornerStorage({std::move(a), std::move(b), std::move(c), std::move(d)}) };
}

template<int dimworld, class Polygon1, class Polygon2>
bool testPolygonIntersection(const Polygon1& pol1,
                             const Polygon2& pol2,
                             bool expectIntersection)
{
    using Policy = Dumux::IntersectionPolicy::PolygonPolicy<double, dimworld>;
    using Algorithm = Dumux::GeometryIntersection<Polygon1, Polygon2, Policy>;

    typename Algorithm::Intersection intersection;
    const bool found = Algorithm::intersection(pol1, pol2, intersection);
    if (found && !expectIntersection)
        std::cout << "Found false positive!" << std::endl;
    else if (!found && expectIntersection)
        std::cout << "Could not find intersection!" << std::endl;
    else if (found && expectIntersection)
        std::cout << "Found intersection" << std::endl;
    else
        std::cout << "No intersection" << std::endl;

    return (found == expectIntersection);
}

template<int dimWorld>
void testPolygonIntersections(std::vector<bool>& returns)
{
    for (auto scaling : {1.0, 1e3, 1e12, 1e-12})
    {
        std::cout << "Test with scaling " << scaling << std::endl;
        const auto tria1 = makeTriangle( makePosition<dimWorld>(0.0, 0.0),
                                         makePosition<dimWorld>(0.0, 1.0*scaling),
                                         makePosition<dimWorld>(1.0*scaling, 1.0*scaling) );
        const auto tria2 = makeTriangle( makePosition<dimWorld>(1.0*scaling, 0.0),
                                         makePosition<dimWorld>(0.0, 0.0),
                                         makePosition<dimWorld>(0.0, 1.0*scaling) );
        const auto tria3 = makeTriangle( makePosition<dimWorld>(1.0*scaling, 0.0),
                                         makePosition<dimWorld>(0.0, 0.5*scaling),
                                         makePosition<dimWorld>(1.0*scaling, 1.0*scaling) );
        const auto tria4 = makeTriangle( makePosition<dimWorld>(1.0*scaling, 0.0),
                                         makePosition<dimWorld>(1.0*scaling, 1.0*scaling),
                                         makePosition<dimWorld>(2.0*scaling, 1.0*scaling) );
        const auto tria5 = makeTriangle( makePosition<dimWorld>(1.0*scaling, 0.0),
                                         makePosition<dimWorld>(0.5*scaling, 0.5*scaling),
                                         makePosition<dimWorld>(2.0*scaling, 1.0*scaling) );
        const auto tria6 = makeTriangle( makePosition<dimWorld>(1.0*scaling, 0.0),
                                         makePosition<dimWorld>(0.501*scaling, 0.501*scaling),
                                         makePosition<dimWorld>(2.0*scaling, 1.0*scaling) );
        const auto tria7 = makeTriangle( makePosition<dimWorld>(1.0*scaling, 0.0),
                                         makePosition<dimWorld>(0.499*scaling, 0.501*scaling),
                                         makePosition<dimWorld>(2.0*scaling, 1.0*scaling) );

        returns.push_back(testPolygonIntersection<dimWorld>(tria1, tria2, true));
        returns.push_back(testPolygonIntersection<dimWorld>(tria1, tria3, true));
        returns.push_back(testPolygonIntersection<dimWorld>(tria1, tria4, false)); // touches in point
        returns.push_back(testPolygonIntersection<dimWorld>(tria1, tria5, false)); // touches in edge
        returns.push_back(testPolygonIntersection<dimWorld>(tria1, tria6, false)); // little bit outside triangle
        returns.push_back(testPolygonIntersection<dimWorld>(tria1, tria7, true));  // little bit inside triangle

        // test some quadrilaterals
        const auto quad1 = makeQuadrilateral( makePosition<dimWorld>(0.0, 0.0),
                                              makePosition<dimWorld>(1.0*scaling, 0.0),
                                              makePosition<dimWorld>(0.0, 1.0*scaling),
                                              makePosition<dimWorld>(1.0*scaling, 1.0*scaling) );
        const auto quad2 = makeQuadrilateral( makePosition<dimWorld>(1.0*scaling, 0.0),
                                              makePosition<dimWorld>(2.0*scaling, 0.0),
                                              makePosition<dimWorld>(1.0*scaling, 2.0*scaling),
                                              makePosition<dimWorld>(2.0*scaling, 2.0*scaling) );
        const auto quad3 = makeQuadrilateral( makePosition<dimWorld>(-1.0*scaling, 0.0),
                                              makePosition<dimWorld>(0.0, 0.0),
                                              makePosition<dimWorld>(-1.0*scaling, 1.0*scaling),
                                              makePosition<dimWorld>(0.0, 1.0*scaling) );
        const auto quad4 = makeQuadrilateral( makePosition<dimWorld>(0.5*scaling, 0.0),
                                              makePosition<dimWorld>(0.5*scaling, 0.501*scaling),
                                              makePosition<dimWorld>(1.0*scaling, 0.0),
                                              makePosition<dimWorld>(1.0*scaling, 1.0*scaling) );

        returns.push_back(testPolygonIntersection<dimWorld>(tria1, quad1, true));
        returns.push_back(testPolygonIntersection<dimWorld>(tria1, quad2, false)); // touches in point
        returns.push_back(testPolygonIntersection<dimWorld>(tria1, quad3, false)); // touches in edge
        returns.push_back(testPolygonIntersection<dimWorld>(tria1, quad4, true)); // small overlap region

        std::cout << std::endl;
    }
}

void testParallelPolygons(std::vector<bool>& returns)
{
    using Point = Dune::FieldVector<double, 3>;

    for (auto scaling : {1.0, 1e3, 1e12, 1e-12})
    {
        const double unit = 1.0*scaling;
        const double offUnit = (1.0 + 1e-6)*scaling;

        std::cout << "Test with scaling " << scaling << std::endl;
        const auto tria1 = makeTriangle( Point{{0.0, 0.0, unit}},
                                         Point{{unit, 0.0, unit}},
                                         Point{{unit, unit, unit}} );
        const auto tria2 = makeTriangle( Point{{0.0, 0.0, offUnit}},
                                         Point{{0.0, unit, offUnit}},
                                         Point{{unit, 0.0, offUnit}} );
        returns.push_back(testPolygonIntersection<3>(tria1, tria2, false));
    }

    std::cout << std::endl;
}

void testNonParallelPolygons(std::vector<bool>& returns)
{
    using Point = Dune::FieldVector<double, 3>;

    for (auto scaling : {1.0, 1e3, 1e12, 1e-12})
    {
        const double unit = 1.0*scaling;
        const double offUnit = (1.0 + 1e-6)*scaling;

        std::cout << "Test with scaling " << scaling << std::endl;
        const auto tria1 = makeTriangle( Point{{0.0, 0.0, unit}},
                                         Point{{unit, 0.0, unit}},
                                         Point{{unit, unit, unit}} );
        const auto tria2 = makeTriangle( Point{{0.0, 0.0, unit}},
                                         Point{{0.0, unit, unit}},
                                         Point{{unit, 0.0, offUnit}} );
        returns.push_back(testPolygonIntersection<3>(tria1, tria2, false));
    }

    std::cout << std::endl;
}

// unit test discovered in a simulation, that originally failed on commit 51fc9c1e3
void testSimulationCase(std::vector<bool>& returns)
{
    using Position = Dune::FieldVector<double, 3>;
    const auto quad1 = makeQuadrilateral(
        Position{{9.49264705847428835739e-01, 1.40073529420942810564e-01, 4.00000000000000077716e-01}},
        Position{{9.44669117627101373458e-01, 1.42095588240539000280e-01, 4.00000000000000077716e-01}},
        Position{{9.49632352923713174420e-01, 1.45036764710471721695e-01, 4.00000000000000022204e-01}},
        Position{{9.46446078418066827354e-01, 1.44730392160359544462e-01, 4.00000000000000077716e-01}}
    );
    const auto quad2 = makeQuadrilateral(
        Position{{9.39999999999999946709e-01, 1.39999999999999985567e-01, 4.00000000000000022204e-01}},
        Position{{9.59999999999999964473e-01, 1.39999999999999985567e-01, 4.00000000000000022204e-01}},
        Position{{9.39999999999999946709e-01, 1.60000000000000003331e-01, 4.00000000000000022204e-01}},
        Position{{9.59999999999999964473e-01, 1.60000000000000003331e-01, 4.00000000000000022204e-01}}
    );

    returns.push_back(testPolygonIntersection<3>(quad1, quad2, true));

    std::cout << std::endl;
}

#endif

int main(int argc, char* argv[])
{
    std::vector<bool> returns;

    std::cout << "Testing intersections in 2d space" << std::endl;
    testPolygonIntersections<2>(returns);

    std::cout << "Testing intersecions in 3d space" << std::endl;
    testPolygonIntersections<3>(returns);

    std::cout << "Testing parallel polygons in 3d space" << std::endl;
    testParallelPolygons(returns);

    std::cout << "Testing non-parallel polygons in 3d space" << std::endl;
    testNonParallelPolygons(returns);

    std::cout << "Testing 3d unit test from actual simulation" << std::endl;
    testSimulationCase(returns);

    // TODO: implement and test point and segment intersections

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), [](bool i){ return !i; }))
        return 1;

    std::cout << "All tests passed!" << std::endl;

    return 0;
}
