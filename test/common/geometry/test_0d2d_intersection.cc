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
#include <dumux/common/geometry/intersectspointgeometry.hh>

#ifndef DOXYGEN
template<class Geometry>
bool testIntersection(const Geometry& geo,
                      const typename Geometry::GlobalCoordinate& p,
                      bool foundExpected, bool verbose)
{
    bool found = Dumux::intersectsPointGeometry(p, geo);
    if (!found && foundExpected)
    {
        std::cerr << "  Failed detecting intersection of " << geo.type();
        for (int i = 0; i < geo.corners(); ++i)
            std::cerr << " (" << geo.corner(i) << ")";
        std::cerr << " with point: " << p << std::endl;
    }
    else if (found && !foundExpected)
    {
        std::cerr << "  Found false positive: intersection of " << geo.type();
        for (int i = 0; i < geo.corners(); ++i)
            std::cerr << " (" << geo.corner(i) << ")";
        std::cerr << " with point: " << p << std::endl;
    }
    if (verbose)
    {
        if (found && foundExpected)
            std::cout << "  Found intersection with " << p << std::endl;
        else if (!found && !foundExpected)
            std::cout << "  No intersection with " << p << std::endl;
    }
    return (found == foundExpected);
}

template<int dimworld, class Transformation>
void runIntersectionTest(std::vector<bool>& returns, const Transformation& transform, bool verbose)
{
    using ctype = typename std::decay_t<decltype(transform({0.0}))>::value_type;
    using Geo = Dune::MultiLinearGeometry<ctype, 2, dimworld>;
    using Points = std::vector<typename Geo::GlobalCoordinate>;

    // test triangle-point intersections
    if (verbose) std::cout << "\n  -- Test triangle-point intersections" << std::endl;

    if constexpr (dimworld == 2)
    {
        auto cornersTri = Points ({{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}});
        std::transform(cornersTri.begin(), cornersTri.end(), cornersTri.begin(),
                       [&transform](const auto& p) { return transform(p); });
        auto triangle = Geo(Dune::GeometryTypes::triangle, cornersTri);

        for (const auto& corner : cornersTri)
            returns.push_back(testIntersection(triangle, corner, true, verbose));

        returns.push_back(testIntersection(triangle, transform({0.25, 0.25}), true, verbose));
        returns.push_back(testIntersection(triangle, transform({0.5, 0.5}), true, verbose));
        returns.push_back(testIntersection(triangle, transform({0.0, 0.5}), true, verbose));
        returns.push_back(testIntersection(triangle, transform({1.01, 0.0}), false, verbose));
        returns.push_back(testIntersection(triangle, transform({0.5, 0.51}), false, verbose));
        returns.push_back(testIntersection(triangle, transform({0.0, -0.01}), false, verbose));
    }

    if constexpr (dimworld == 3)
    {
        auto cornersTri = Points ({{0.0, 0.0, 2.0}, {1.0, 0.0, 2.0}, {0.0, 1.0, 2.0}});
        std::transform(cornersTri.begin(), cornersTri.end(), cornersTri.begin(),
                       [&transform](const auto& p) { return transform(p); });
        auto triangle = Geo(Dune::GeometryTypes::triangle, cornersTri);

        for (const auto& corner : cornersTri)
            returns.push_back(testIntersection(triangle, corner, true, verbose));

        returns.push_back(testIntersection(triangle, transform({0.25, 0.25, 2.0}), true, verbose));
        returns.push_back(testIntersection(triangle, transform({0.5, 0.5, 2.0}), true, verbose));
        returns.push_back(testIntersection(triangle, transform({0.0, 0.5, 2.0}), true, verbose));
        returns.push_back(testIntersection(triangle, transform({1.01, 0.0, 2.0}), false, verbose));
        returns.push_back(testIntersection(triangle, transform({0.5, 0.51, 2.0}), false, verbose));
        returns.push_back(testIntersection(triangle, transform({0.0, -0.01, 2.0}), false, verbose));
    }

    // test quadrilateral-point intersections
    if (verbose) std::cout << "\n  -- Test quadrilateral-point intersections" << std::endl;

    if constexpr (dimworld == 2)
    {
        auto cornersQuad = Points ({{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}});
        std::transform(cornersQuad.begin(), cornersQuad.end(), cornersQuad.begin(),
                       [&transform](const auto& p) { return transform(p); });
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

    if constexpr (dimworld == 3)
    {
        auto cornersQuad = Points ({{0.0, 0.0, 2.0}, {1.0, 0.0, 2.0}, {0.0, 1.0, 2.0}, {1.0, 1.0, 2.0}});
        std::transform(cornersQuad.begin(), cornersQuad.end(), cornersQuad.begin(),
                       [&transform](const auto& p) { return transform(p); });
        auto quadrilateral = Geo(Dune::GeometryTypes::quadrilateral, cornersQuad);

        for (const auto& corner : cornersQuad)
            returns.push_back(testIntersection(quadrilateral, corner, true, verbose));

        returns.push_back(testIntersection(quadrilateral, transform({0.5, 0.5, 2.0}), true, verbose));
        returns.push_back(testIntersection(quadrilateral, transform({0.5, 0.0, 2.0}), true, verbose));
        returns.push_back(testIntersection(quadrilateral, transform({0.5, 1.0, 2.0}), true, verbose));
        returns.push_back(testIntersection(quadrilateral, transform({1.01, 1.0, 2.0}), false, verbose));
        returns.push_back(testIntersection(quadrilateral, transform({0.5, 1.01, 2.0}), false, verbose));
        returns.push_back(testIntersection(quadrilateral, transform({0.0, -0.01, 2.0}), false, verbose));
    }
}

template<class ctype>
auto create2DTransformation(const ctype scale,
                          const Dune::FieldVector<ctype, 2>& translate,
                          const ctype rotationAngle)
{
    std::cout << "Intersection test with transformation:"
              << " ctype: " << Dune::className<ctype>()
              << ", scaling: " << scale
              << ", translation: " << translate
              << ", rotationAngle: " << rotationAngle << std::endl;
    // Rotation of a vector in two dimensions
    // See rotation matrix (2d) at https://en.wikipedia.org/wiki/Rotation_matrix
    const ctype sinAngle = std::sin(rotationAngle);
    const ctype cosAngle = std::cos(rotationAngle);
    return [=](Dune::FieldVector<ctype, 2> p){
        p *= scale;
        p.axpy(scale, translate);
        auto tp = p;
        tp[0] = p[0]*cosAngle-p[1]*sinAngle;
        tp[1] = p[0]*sinAngle+p[1]*cosAngle;
        return tp;
    };
}

template<class ctype>
auto create3DTransformation(const ctype scale,
                          const Dune::FieldVector<ctype, 3>& translate,
                          const Dune::FieldVector<ctype, 3>& rotationAxis,
                          const ctype rotationAngle)
{
    std::cout << "Intersection test with transformation:"
              << " ctype: " << Dune::className<ctype>()
              << ", scaling: " << scale
              << ", translation: " << translate
              << ", rotationAxis: " << rotationAxis
              << ", rotationAngle: " << rotationAngle << std::endl;
    // Rotation of a vector in three dimensions
    // See Rodrigues' rotation formular at https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    const ctype sinAngle = std::sin(rotationAngle);
    const ctype cosAngle = std::cos(rotationAngle);
    return [=](Dune::FieldVector<ctype, 3> p){
        p *= scale;
        p.axpy(scale, translate);
        auto tp = p;
        tp *= cosAngle;
        tp.axpy(sinAngle, Dumux::crossProduct({rotationAxis}, p));
        return tp.axpy((1.0-cosAngle)*(rotationAxis*p), rotationAxis);
    };
}

#endif

int main (int argc, char *argv[]) try
{
    // collect returns to determine exit code
    std::vector<bool> returns;
    constexpr bool verbose = false;

    {
        constexpr int dimworld = 2;
        using Vec = Dune::FieldVector<double, dimworld>;

        std::cout << "Test for point intersection with 2d geometries in 2d corddim" << std::endl;
        for (const double scaling : {1.0, 1e3, 1e12, 1e-12})
            for (const auto& translation : {0.0, 1.0})
                for (const double angle : {0.0, 0.2*M_PI, 0.5*M_PI, 0.567576567*M_PI, M_PI})
                    runIntersectionTest<dimworld>(returns, create2DTransformation<double>(scaling, Vec(translation), angle), verbose);
    }

    {
        constexpr int dimworld = 3;
        using Vec = Dune::FieldVector<double, dimworld>;

        std::cout << "Test for point intersection with 2d geometries in 3d corddim" << std::endl;
        for (const double scaling : {1.0, 1e3, 1e12, 1e-12})
            for (const auto& translation : {0.0, 1.0})
                for (const double angle : {0.0, 0.2*M_PI, 0.5*M_PI, 0.567576567*M_PI, M_PI})
                    for (const auto& rotAxis : {Vec(std::sqrt(3.0)/3.0), Vec({std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0})})
                        runIntersectionTest<dimworld>(returns, create3DTransformation<double>(scaling, Vec(translation), rotAxis, angle), verbose);
    }

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), std::logical_not<bool>{}))
        return 1;

    std::cout << "\n++++++++++++++++++++++\n"
              << "All tests passed!"
              << "\n++++++++++++++++++++++" << std::endl;

    return 0;
}
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
