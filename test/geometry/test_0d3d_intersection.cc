#include <config.h>

#include <iostream>
#include <algorithm>
#include <functional>
#include <type_traits>

#include <dune/common/classname.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dumux/common/math.hh>

#include "transformation.hh"
#include "test_intersection.hh"

namespace Dumux {

template<class Transformation>
void runIntersectionTest(std::vector<bool>& returns, const Transformation& transform, bool verbose)
{
    using ctype = typename std::decay_t<decltype(transform({0.0}))>::value_type;
    using Geo = Dune::MultiLinearGeometry<ctype, 3, 3>;
    using Points = std::vector<typename Geo::GlobalCoordinate>;

    // test tetrahedron-point intersections
    if (verbose) std::cout << "\n  -- Test tetrahedron-point intersections" << std::endl;

    auto cornersTet = Points({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
    std::transform(cornersTet.begin(), cornersTet.end(), cornersTet.begin(),
                   [&](const auto& p) { return transform(p); });
    const auto tetrahedron = Geo(Dune::GeometryTypes::tetrahedron, cornersTet);

    for (const auto& corner : cornersTet)
        returns.push_back(testIntersection(tetrahedron, corner, true, verbose));
    returns.push_back(testIntersection(tetrahedron, transform({0.0, 0.0, 0.5}), true, verbose));
    returns.push_back(testIntersection(tetrahedron, transform({0.25, 0.25, 0.5}), true, verbose));
    returns.push_back(testIntersection(tetrahedron, transform({0.5, 0.5, 0.5}), false, verbose));
    returns.push_back(testIntersection(tetrahedron, transform({1.01, 0.0, 0.0}), false, verbose));
    returns.push_back(testIntersection(tetrahedron, transform({0.5, 0.0, 0.51}), false, verbose));

    // test hexahedron-point intersections
    if (verbose) std::cout << "\n  -- Test hexahedron-point intersections" << std::endl;

    auto cornersHex = Points({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {1.0, 1.0, 0.0},
                              {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}});
    std::transform(cornersHex.begin(), cornersHex.end(), cornersHex.begin(),
                   [&](const auto& p) { return transform(p); });
    auto hexahedron = Geo(Dune::GeometryTypes::hexahedron, cornersHex);

    for (const auto& corner : cornersHex)
        returns.push_back(testIntersection(hexahedron, corner, true, verbose));
    returns.push_back(testIntersection(hexahedron, transform({0.5, 0.5, 0.5}), true, verbose));
    returns.push_back(testIntersection(hexahedron, transform({0.01, 0.01, 0.0001}), true, verbose));
    returns.push_back(testIntersection(hexahedron, transform({1.01, 0.5, 0.5}), false, verbose));
    returns.push_back(testIntersection(hexahedron, transform({2.0, 2.0, 2.0}), false, verbose));
    returns.push_back(testIntersection(hexahedron, transform({-0.5, -0.0, -0.51}), false, verbose));

    // test pyramid-point intersections
    if (verbose) std::cout << "\n  -- Test pyramid-point intersections" << std::endl;

    auto cornersPyramid = Points({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                                  {1.0, 1.0, 0.0}, {0.5, 0.5, 1.0}});
    std::transform(cornersPyramid.begin(), cornersPyramid.end(), cornersPyramid.begin(),
                   [&](const auto& p) { return transform(p); });
    auto pyramid = Geo(Dune::GeometryTypes::pyramid, cornersPyramid);

    for (const auto& corner : cornersPyramid)
        returns.push_back(testIntersection(pyramid, corner, true, verbose));
    returns.push_back(testIntersection(pyramid, transform({0.5, 0.5, 0.0}), true, verbose));
    returns.push_back(testIntersection(pyramid, transform({0.5, 0.5, 0.7}), true, verbose));
    returns.push_back(testIntersection(pyramid, transform({0.5, 0.5, -0.0001}), false, verbose));
    returns.push_back(testIntersection(pyramid, transform({0.25, 0.25, 0.5}), true, verbose));
    returns.push_back(testIntersection(pyramid, transform({0.25, 0.75, 0.5}), true, verbose));
    returns.push_back(testIntersection(pyramid, transform({0.75, 0.75, 0.5}), true, verbose));
    returns.push_back(testIntersection(pyramid, transform({0.25, 0.25, 0.5001}), false, verbose));
    returns.push_back(testIntersection(pyramid, transform({1.0, 1.0, 0.0001}), false, verbose));

    // test prism-point intersections
    if (verbose) std::cout << "\n  -- Test prism-point intersections" << std::endl;

    auto cornersPrism = Points({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                                {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}});
    std::transform(cornersPrism.begin(), cornersPrism.end(), cornersPrism.begin(),
                   [&](const auto& p) { return transform(p); });
    auto prism = Geo(Dune::GeometryTypes::prism, cornersPrism);

    for (const auto& corner : cornersPrism)
        returns.push_back(testIntersection(prism, corner, true, verbose));
    returns.push_back(testIntersection(prism, transform({0.25, 0.0, 0.25}), true, verbose));
    returns.push_back(testIntersection(prism, transform({0.0, 0.25, 0.25}), true, verbose));
    returns.push_back(testIntersection(prism, transform({0.25, 0.25, 1.0}), true, verbose));
    returns.push_back(testIntersection(prism, transform({0.25, 0.25, 0.0}), true, verbose));
    returns.push_back(testIntersection(prism, transform({0.25, 0.25, -0.0001}), false, verbose));
    returns.push_back(testIntersection(prism, transform({0.25, 0.25, 1.0001}), false, verbose));
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
                    runIntersectionTest(returns,
                        make3DTransformation<double>(scaling, Vec(translation), rotAxis, angle, true), verbose);

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), std::logical_not<bool>{}))
        return 1;

    std::cout << "\n++++++++++++++++++++++\n"
              << "All tests passed!"
              << "\n++++++++++++++++++++++" << std::endl;

    return 0;
}
