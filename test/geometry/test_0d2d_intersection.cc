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

template<int dimworld, class Transformation>
void runIntersectionTest(std::vector<bool>& returns, const Transformation& transform, bool verbose)
{
    using ctype = typename std::decay_t<decltype(transform({0.0}))>::value_type;
    using Geo = Dune::MultiLinearGeometry<ctype, 2, dimworld>;
    using Points2D = std::vector<Dune::FieldVector<ctype, 2>>;
    using PointsDimWorld = std::vector<Dune::FieldVector<ctype, dimworld>>;

    // test triangle-point intersections
    if (verbose) std::cout << "\n  -- Test triangle-point intersections" << std::endl;

    auto cornersTri2D = Points2D ({{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}});
    auto cornersTri = PointsDimWorld(3);
    std::transform(cornersTri2D.begin(), cornersTri2D.end(), cornersTri.begin(),
                   [&](const auto& p) { return transform(p); });
    auto triangle = Geo(Dune::GeometryTypes::triangle, cornersTri);

    for (const auto& corner : cornersTri)
        returns.push_back(testIntersection(triangle, corner, true, verbose));

    returns.push_back(testIntersection(triangle, transform({0.25, 0.25}), true, verbose));
    returns.push_back(testIntersection(triangle, transform({0.5, 0.5}), true, verbose));
    returns.push_back(testIntersection(triangle, transform({0.0, 0.5}), true, verbose));
    returns.push_back(testIntersection(triangle, transform({1.01, 0.0}), false, verbose));
    returns.push_back(testIntersection(triangle, transform({0.5, 0.51}), false, verbose));
    returns.push_back(testIntersection(triangle, transform({0.0, -0.01}), false, verbose));

    // test quadrilateral-point intersections
    if (verbose) std::cout << "\n  -- Test quadrilateral-point intersections" << std::endl;

    auto cornersQuad2D = Points2D ({{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}});
    auto cornersQuad = PointsDimWorld(4);
    std::transform(cornersQuad2D.begin(), cornersQuad2D.end(), cornersQuad.begin(),
                   [&](const auto& p) { return transform(p); });
    auto quadrilateral = Geo(Dune::GeometryTypes::quadrilateral, cornersQuad);

    for (const auto& corner : cornersQuad)
        returns.push_back(testIntersection(quadrilateral, corner, true, verbose));

    returns.push_back(testIntersection(quadrilateral, transform({0.5, 0.5}), true, verbose));
    returns.push_back(testIntersection(quadrilateral, transform({0.5, 0.0}), true, verbose));
    returns.push_back(testIntersection(quadrilateral, transform({0.5, 1.0}), true, verbose));
    returns.push_back(testIntersection(quadrilateral, transform({1.01, 1.0}), false, verbose));
    returns.push_back(testIntersection(quadrilateral, transform({0.5, 1.01}), false, verbose));
    returns.push_back(testIntersection(quadrilateral, transform({0.0, -0.01}), false, verbose));

}

template<class Transformation>
void run3DIntersectionTest(std::vector<bool>& returns, const Transformation& transform, bool verbose)
{
    // augment 2d vector by a third coordinate
    using Point = std::decay_t<decltype(transform({0.0, 0.0, 0.0}))>;
    const auto transform2D3D = [&](Dune::FieldVector<typename Point::value_type, 2> p2D){
        return transform({p2D[0], p2D[1], 2.0});
    };

    runIntersectionTest<3>(returns, transform2D3D, verbose);
}

template<class Transformation>
void run2DIntersectionTest(std::vector<bool>& returns, const Transformation& transform, bool verbose)
{ runIntersectionTest<2>(returns, transform, verbose); }

} // end namespace Dumux

int main (int argc, char *argv[])
{
    using namespace Dumux;

    // collect returns to determine exit code
    std::vector<bool> returns;
    constexpr bool verbose = false;

    {
        using Vec = Dune::FieldVector<double, 2>;

        std::cout << "Test for point intersection with 2d geometries in 2d corddim" << std::endl;
        for (const double scaling : {1.0, 1e3, 1e12, 1e-12})
            for (const auto& translation : {0.0, 1.0})
                for (const double angle : {0.0, 0.2*M_PI, 0.5*M_PI, 0.567576567*M_PI, M_PI})
                    run2DIntersectionTest(returns,
                        make2DTransformation<double>(scaling, Vec(translation), angle, true), verbose);
    }

    {
        using Vec = Dune::FieldVector<double, 3>;

        std::cout << "Test for point intersection with 2d geometries in 3d corddim" << std::endl;
        for (const double scaling : {1.0, 1e3, 1e12, 1e-12})
            for (const auto& translation : {0.0, 1.0})
                for (const double angle : {0.0, 0.2*M_PI, 0.5*M_PI, 0.567576567*M_PI, M_PI})
                    for (const auto& rotAxis : {Vec(std::sqrt(3.0)/3.0), Vec({std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0})})
                        run3DIntersectionTest(returns,
                            make3DTransformation<double>(scaling, Vec(translation), rotAxis, angle, true), verbose);
    }

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), std::logical_not<bool>{}))
        return 1;

    std::cout << "\n++++++++++++++++++++++\n"
              << "All tests passed!"
              << "\n++++++++++++++++++++++" << std::endl;

    return 0;
}
