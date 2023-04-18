//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>

#include "transformation.hh"
#include <dumux/geometry/geometryintersection.hh>

#ifndef DOXYGEN
Dune::MultiLinearGeometry<double, 1, 3>
makeLine(std::initializer_list<Dune::FieldVector<double, 3>>&& c)
{
    return {Dune::GeometryTypes::line, c};
}

bool testIntersection(const Dune::MultiLinearGeometry<double, 3, 3>& polyhedron,
                      const Dune::MultiLinearGeometry<double, 1, 3>& line,
                      bool foundExpected, bool verbose)
{
    using Test = Dumux::GeometryIntersection<Dune::MultiLinearGeometry<double, 3, 3>,
                                             Dune::MultiLinearGeometry<double, 1, 3> >;
    typename Test::Intersection intersection;
    bool found = Test::intersection(polyhedron, line, intersection);
    if (!found && foundExpected)
    {
        std::cerr << "  Failed detecting intersection with " << line.corner(0) << " " << line.corner(1) << std::endl;
    }
    else if (found && !foundExpected)
    {
        std::cerr << "  Found false positive: intersection with " << line.corner(0) << " " << line.corner(1) << std::endl;
    }
    if (verbose)
    {
        if (found && foundExpected)
            std::cout << "  Found intersection with " << line.corner(0) << " " << line.corner(1) << std::endl;
        else if (!found && !foundExpected)
            std::cout << "  No intersection with " << line.corner(0) << " " << line.corner(1) << std::endl;
    }
    return (found == foundExpected);
}
#endif

namespace Dumux {

void testEdgesAndDiagonals(std::vector<bool>& returns,
                           const Dune::MultiLinearGeometry<double, 3, 3>& polyhedron,
                           bool verbose)
{
    for (int i = 0; i < polyhedron.corners(); ++i)
    {
        for (int j = 0; j < polyhedron.corners(); ++j)
        {
            const auto& ci = polyhedron.corner(i);
            const auto& cj = polyhedron.corner(j);
            returns.push_back(testIntersection(polyhedron, makeLine({ci, cj}), true, verbose));
            returns.push_back(testIntersection(polyhedron, makeLine({ci, cj + (ci-cj)}), true, verbose));
            returns.push_back(testIntersection(polyhedron, makeLine({ci - (ci-cj), cj}), true, verbose));
        }
    }
}

template<class Transformation>
void runIntersectionTest(std::vector<bool>& returns, const Transformation& transform, bool verbose)
{
    using Points = std::vector<Dune::FieldVector<double, 3>>;
    using Geo = Dune::MultiLinearGeometry<double, 3, 3>;

    // test tetrahedron-line intersections
    if (verbose) std::cout << "\n  -- Test tetrahedron-line intersections" << std::endl;

    auto cornersTetrahedron = Points({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
    std::transform(cornersTetrahedron.begin(), cornersTetrahedron.end(), cornersTetrahedron.begin(),
                   [&](const auto& p) { return transform(p); });
    const auto tetrahedron = Geo(Dune::GeometryTypes::tetrahedron, cornersTetrahedron);

    // test all edges and diagonals in both orientations
    testEdgesAndDiagonals(returns, tetrahedron, verbose);

    returns.push_back(testIntersection(tetrahedron, makeLine({{transform({0.25, 0.25, 0.0})}, {transform({0.25, 0.25, 1.0})}}), true, verbose));
    returns.push_back(testIntersection(tetrahedron, makeLine({{transform({-1.0, 0.25, 0.5})}, {transform({1.0, 0.25, 0.5})}}), true, verbose));
    returns.push_back(testIntersection(tetrahedron, makeLine({{transform({1.0, 1.0, 1.0})}, {transform({-1.0, -1.0, -1.0})}}), true, verbose));

    returns.push_back(testIntersection(tetrahedron, makeLine({{transform({1.5, 0.0, 0.5})}, {transform({0.0, 1.5, 0.5})}}), false, verbose));
    returns.push_back(testIntersection(tetrahedron, makeLine({{transform({0.0, 0.0, 0.0})}, {transform({0.0, 0.0, -1.0})}}), false, verbose));
    returns.push_back(testIntersection(tetrahedron, makeLine({{transform({1.0, 1.0, 0.0})}, {transform({0.0, 0.0, 2.0})}}), false, verbose));
    returns.push_back(testIntersection(tetrahedron, makeLine({{transform({1.0, 0.0, 0.1})}, {transform({0.0, 1.0, 0.1})}}), false, verbose));
    returns.push_back(testIntersection(tetrahedron, makeLine({{transform({0.0, 0.0, -0.1})}, {transform({1.0, 1.0, -0.1})}}), false, verbose));

    // lines with a small offset
    returns.push_back(testIntersection(tetrahedron, makeLine({{transform({0.0, -1e-5, 0.0})}, {transform({1.0, 0.0, 0.0})}}), false, verbose));
    returns.push_back(testIntersection(tetrahedron, makeLine({{transform({0.0, -1e-5, 0.0})}, {transform({1.0, -1e-5, 0.0})}}), false, verbose));
    returns.push_back(testIntersection(tetrahedron, makeLine({{transform({0.0, 0.0, 0.0})}, {transform({1.0, -1e-5, 0.0})}}), false, verbose));

    // test hexahedron-line intersections
    if (verbose) std::cout << "\n  -- Test hexahedron-line intersections" << std::endl;

    auto cornersHexahedron = Points({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0},
                                     {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}});
    std::transform(cornersHexahedron.begin(), cornersHexahedron.end(), cornersHexahedron.begin(),
                   [&](const auto& p) { return transform(p); });
    auto hexahedron = Geo(Dune::GeometryTypes::hexahedron, cornersHexahedron);

    // test all edges and diagonals in both orientations
    testEdgesAndDiagonals(returns, hexahedron, verbose);

    returns.push_back(testIntersection(hexahedron, makeLine({{transform({0.0, 0.0, 2.0})}, {transform({1.0, 1.0, 2.0})}}), false, verbose));
    returns.push_back(testIntersection(hexahedron, makeLine({{transform({0.0, 0.0, 1.1})}, {transform({1.0, 1.0, 1.1})}}), false, verbose));
    returns.push_back(testIntersection(hexahedron, makeLine({{transform({1.1, 1.1, 0.0})}, {transform({1.1, 1.1, 1.0})}}), false, verbose));
    returns.push_back(testIntersection(hexahedron, makeLine({{transform({1.1, 0.0, 0.0})}, {transform({1.1, 1.0, 1.0})}}), false, verbose));
    returns.push_back(testIntersection(hexahedron, makeLine({{transform({0.0, -0.1, 0.0})}, {transform({1.0, -0.1, 0.0})}}), false, verbose));

    // test pyramid-line intersections
    if (verbose) std::cout << "\n  -- Test pyramid-line intersections" << std::endl;

    auto cornersPyramid = Points({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0}, {0.5, 0.5, 1.0}});
    std::transform(cornersPyramid.begin(), cornersPyramid.end(), cornersPyramid.begin(),
                   [&](const auto& p) { return transform(p); });
    auto pyramid = Geo(Dune::GeometryTypes::pyramid, cornersPyramid);

    // test all edges and diagonals in both orientations
    testEdgesAndDiagonals(returns, pyramid, verbose);

    returns.push_back(testIntersection(pyramid, makeLine({{transform({0.0, 0.0, 1.0})}, {transform({1.0, 0.0, 1.0})}}), false, verbose));
    returns.push_back(testIntersection(pyramid, makeLine({{transform({0.0, 0.0, 1.0})}, {transform({0.0, 1.0, 1.0})}}), false, verbose));
    returns.push_back(testIntersection(pyramid, makeLine({{transform({0.0, 0.0, -0.1})}, {transform({1.0, 1.0, -0.1})}}), false, verbose));
    returns.push_back(testIntersection(pyramid, makeLine({{transform({0.0, 1.1, 0.0})}, {transform({1.0, 1.1, 0.0})}}), false, verbose));
    returns.push_back(testIntersection(pyramid, makeLine({{transform({0.4, 0.0, 1.0})}, {transform({0.4, 1.0, 1.0})}}), false, verbose));

    // test prism-line intersections
    if (verbose) std::cout << "\n  -- Test prism-line intersections" << std::endl;

    auto cornersPrism = Points({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                                {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}});
    std::transform(cornersPrism.begin(), cornersPrism.end(), cornersPrism.begin(),
                   [&](const auto& p) { return transform(p); });
    auto prism = Geo(Dune::GeometryTypes::prism, cornersPrism);

    // test all edges and diagonals in both orientations
    testEdgesAndDiagonals(returns, prism, verbose);

    returns.push_back(testIntersection(prism, makeLine({{transform({1.0, 1.0, 0.0})}, {transform({1.0, 1.0, 1.0})}}), false, verbose));
    returns.push_back(testIntersection(prism, makeLine({{transform({2.0, 0.0, 0.5})}, {transform({0.0, 2.0, 0.5})}}), false, verbose));
    returns.push_back(testIntersection(prism, makeLine({{transform({1.1, 0.0, 0.0})}, {transform({1.1, 0.0, 1.0})}}), false, verbose));
    returns.push_back(testIntersection(prism, makeLine({{transform({-0.1, 0.0, 1.0})}, {transform({-0.1, 1.0, 1.0})}}), false, verbose));
    returns.push_back(testIntersection(prism, makeLine({{transform({1.0, 0.0, 1.1})}, {transform({0.0, 1.0, 1.1})}}), false, verbose));
}

} // end namespace Dumux

int main (int argc, char *argv[])
{
    using namespace Dumux;

    // collect returns to determine exit code
    std::vector<bool> returns;
    constexpr bool verbose = false;

    using Vec = Dune::FieldVector<double, 3>;
    for (const double scaling : {1.0, 1e3, 1e12, 1e-12})
        for (const double translation : {0.0, 1.0})
            for (const double angle : {0.0, 0.2*M_PI, 0.5*M_PI, 0.567576567*M_PI, M_PI})
                for (const auto& rotAxis : {Vec(std::sqrt(3.0)/3.0), Vec({std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0})})
                    runIntersectionTest(returns, make3DTransformation<double>(scaling, Vec(translation), rotAxis, angle, true), verbose);

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), std::logical_not<bool>{}))
        return 1;

    std::cout << "\n++++++++++++++++++++++\n"
              << "All tests passed!"
              << "\n++++++++++++++++++++++" << std::endl;

    return 0;
}
