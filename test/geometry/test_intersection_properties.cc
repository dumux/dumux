//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Geometry
 * \brief Properties of the geometry intersections for random convex geometries
 *
 * Random convex geometries of all supported element types are intersected with the cells of
 * structured grids, with each other and with segments, under several rotations, scalings and
 * translations. Checked are: the measures of the intersections with all cells of a grid add up to
 * the measure of the geometry, both argument orders give the same intersection, point-like
 * intersections agree with an independent computation, tree queries agree with testing all pairs,
 * and cells of a conforming grid only intersect if they are moved into each other.
 */
#include <config.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/math.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/geometryintersection.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/geometry/triangulation.hh>

namespace {

template<int dim> using Point = Dune::FieldVector<double, dim>;
template<int mydim, int dim> using Geo = Dune::MultiLinearGeometry<double, mydim, dim>;

std::mt19937_64 rng(1);
double uniform(double a, double b) { return a + (b - a)*std::generate_canonical<double, 64>(rng); }

////////////////////////////////////////////////////////////////////////////////
// results
////////////////////////////////////////////////////////////////////////////////

int numFailures = 0;

void report(const std::string& what, double error, double tolerance, std::size_t samples)
{
    const bool ok = error <= tolerance;
    if (!ok)
        ++numFailures;
    std::cout << (ok ? "  ok   " : "  FAIL ") << what << ": max error " << error
              << " (tolerance " << tolerance << ", " << samples << " samples)" << std::endl;
}

void reportCount(const std::string& what, std::size_t wrong, std::size_t samples)
{
    if (wrong > 0)
        ++numFailures;
    std::cout << (wrong == 0 ? "  ok   " : "  FAIL ") << what << ": " << wrong << " wrong of "
              << samples << " samples" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
// measures
////////////////////////////////////////////////////////////////////////////////

template<int dim>
Point<3> toPoint3(const Point<dim>& p)
{
    Point<3> q(0.0);
    for (int i = 0; i < dim; ++i)
        q[i] = p[i];
    return q;
}

template<int dim>
double simplexMeasure(const Point<dim>& a, const Point<dim>& b)
{ return (b - a).two_norm(); }

template<int dim>
double simplexMeasure(const Point<dim>& a, const Point<dim>& b, const Point<dim>& c)
{ return 0.5*Dumux::crossProduct(toPoint3(b) - toPoint3(a), toPoint3(c) - toPoint3(a)).two_norm(); }

double simplexMeasure(const Point<3>& a, const Point<3>& b, const Point<3>& c, const Point<3>& d)
{ using std::abs; return abs(Dumux::tripleProduct(b - a, c - a, d - a))/6.0; }

//! measure of a simplex given by a random-access container of corners
template<class Corners>
double simplexMeasure(const Corners& c)
{
    if (c.size() == 1) return 0.0;
    if (c.size() == 2) return simplexMeasure(c[0], c[1]);
    if (c.size() == 3) return simplexMeasure(c[0], c[1], c[2]);
    if constexpr (std::decay_t<decltype(c[0])>::dimension == 3)
        if (c.size() == 4) return simplexMeasure(c[0], c[1], c[2], c[3]);
    DUNE_THROW(Dune::InvalidStateException, "No simplex with " << c.size() << " corners");
}

//! measure of a convex element as the sum of a subdivision into simplices (pulling from corner 0)
template<int mydim, int dim>
double elementMeasure(const Geo<mydim, dim>& g)
{
    const auto c = [&](int i){ return g.corner(i); };
    const auto type = g.type();
    if (type.isLine()) return simplexMeasure(c(0), c(1));
    if (type.isTriangle()) return simplexMeasure(c(0), c(1), c(2));
    if (type.isQuadrilateral()) return simplexMeasure(c(0), c(1), c(3)) + simplexMeasure(c(0), c(3), c(2));
    if constexpr (dim == 3)
    {
        if (type.isTetrahedron()) return simplexMeasure(c(0), c(1), c(2), c(3));
        if (type.isPyramid()) return simplexMeasure(c(0), c(1), c(3), c(4)) + simplexMeasure(c(0), c(3), c(2), c(4));
        if (type.isPrism()) return simplexMeasure(c(0), c(1), c(2), c(5)) + simplexMeasure(c(0), c(1), c(5), c(4))
                                   + simplexMeasure(c(0), c(3), c(4), c(5));
        if (type.isHexahedron()) return simplexMeasure(c(0), c(4), c(5), c(7)) + simplexMeasure(c(0), c(4), c(7), c(6))
                                        + simplexMeasure(c(0), c(2), c(3), c(7)) + simplexMeasure(c(0), c(2), c(7), c(6))
                                        + simplexMeasure(c(0), c(1), c(3), c(7)) + simplexMeasure(c(0), c(1), c(7), c(5));
    }
    DUNE_THROW(Dune::InvalidStateException, "Unknown element type " << type);
}

//! measure of the result of GeometryIntersection
template<class Policy, class Intersection>
double intersectionMeasure(const Intersection& is)
{
    constexpr int dimIs = Policy::dimIntersection;
    constexpr int dimWorld = Policy::dimWorld;
    if constexpr (dimIs == 0) return 0.0;
    else if constexpr (dimIs == 1) return simplexMeasure(is[0], is[1]);
    else
    {
        double m = 0.0;
        for (const auto& s : Dumux::triangulate<dimIs, dimWorld>(is))
            m += simplexMeasure(s);
        return m;
    }
}

template<class Container>
double sumOfPieces(const Container& intersections)
{
    double m = 0.0;
    for (const auto& is : intersections)
        m += simplexMeasure(is.corners());
    return m;
}

////////////////////////////////////////////////////////////////////////////////
// transformations
////////////////////////////////////////////////////////////////////////////////

template<int dim>
std::array<Point<dim>, dim> randomRotation()
{
    if constexpr (dim == 2)
    {
        const double a = uniform(0.0, 2.0*M_PI);
        return {{{std::cos(a), -std::sin(a)}, {std::sin(a), std::cos(a)}}};
    }
    else
    {
        double q[4]; double n = 0.0;
        for (auto& x : q) { x = uniform(-1.0, 1.0); n += x*x; }
        n = std::sqrt(n);
        for (auto& x : q) x /= n;
        const double w = q[0], x = q[1], y = q[2], z = q[3];
        return {{{1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)},
                 {2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)},
                 {2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)}}};
    }
}

template<class Rows, int dim>
Point<dim> rotate(const Rows& R, const Point<dim>& v)
{
    Point<dim> r;
    for (int i = 0; i < dim; ++i)
        r[i] = R[i]*v;
    return r;
}

//! x -> scale*(R x + shift), the shift scales with the geometry like in a change of units,
//! R is a rotation or a rotation combined with a reflection, which reverses the orientation of all elements
template<int dim>
struct Transformation
{
    std::array<Point<dim>, dim> rotation;
    double scale;
    Point<dim> shift;
    bool mirrored = false;

    Point<dim> operator()(const Point<dim>& x) const
    {
        auto y = rotate(rotation, x);
        y += shift;
        y *= scale;
        return y;
    }

    std::string name() const
    {
        std::ostringstream s;
        s << "scale " << scale << ", shift " << shift.infinity_norm() << (mirrored ? ", mirrored" : "");
        return s.str();
    }
};

template<int dim>
std::vector<Transformation<dim>> transformations()
{
    std::vector<Transformation<dim>> t;
    for (const auto& [scale, distance] : {std::pair{1.0, 2.0}, {1e-5, 2.0}, {1e5, 2.0}, {1.0, 1e3}})
    {
        for (const bool mirrored : {false, true})
        {
            Point<dim> shift;
            for (auto& s : shift)
                s = uniform(-distance, distance);
            auto rotation = randomRotation<dim>();
            if (mirrored)
                rotation[0] *= -1.0;
            t.push_back({rotation, scale, shift, mirrored});
        }
    }
    return t;
}

template<int mydim, int dim, class Map>
Geo<mydim, dim> makeGeometry(Dune::GeometryType type, const std::vector<Point<dim>>& corners, const Map& map)
{
    std::vector<Point<dim>> c;
    for (const auto& p : corners)
        c.push_back(map(p));
    return Geo<mydim, dim>(type, c);
}

//! geometry in 3d from corners in a plane and the map of the plane into 3d
template<int mydim, class Map>
Geo<mydim, 3> makeEmbeddedGeometry(Dune::GeometryType type, const std::vector<Point<2>>& corners, const Map& map)
{
    std::vector<Point<3>> c;
    for (const auto& p : corners)
        c.push_back(map(p));
    return Geo<mydim, 3>(type, c);
}

////////////////////////////////////////////////////////////////////////////////
// random convex elements (reference coordinates, before the transformation)
////////////////////////////////////////////////////////////////////////////////

//! points on an ellipse in cyclic order, around the origin
std::vector<Point<2>> ellipsePoints(int n)
{
    const double rx = uniform(0.6, 1.0), ry = uniform(0.6, 1.0), phase = uniform(0.0, 2.0*M_PI);
    std::vector<Point<2>> p;
    for (int k = 0; k < n; ++k)
    {
        const double phi = phase + 2.0*M_PI*(k + uniform(-0.25, 0.25))/n;
        p.push_back({rx*std::cos(phi), ry*std::sin(phi)});
    }
    return p;
}

//! corners in Dune order of a random convex polygon (triangle or quadrilateral) of size r around c
template<int dim>
std::vector<Point<dim>> randomPolygon(Dune::GeometryType type, const Point<dim>& c, double r,
                                      const std::array<Point<dim>, dim>& R)
{
    auto p = ellipsePoints(type.isTriangle() ? 3 : 4);
    if (type.isQuadrilateral())
        std::swap(p[2], p[3]);
    std::vector<Point<dim>> corners;
    for (const auto& q : p)
    {
        Point<dim> x(0.0);
        x[0] = r*q[0]; x[1] = r*q[1];
        corners.push_back(c + rotate(R, x));
    }
    return corners;
}

//! corners in Dune order of a random convex polyhedron with planar faces of size r around c
std::vector<Point<3>> randomPolyhedron(Dune::GeometryType type, const Point<3>& c, double r)
{
    const auto R = randomRotation<3>();
    const auto place = [&](const Point<3>& x){ return c + rotate(R, r*x); };
    std::vector<Point<3>> corners;
    if (type.isTetrahedron())
    {
        while (true)
        {
            corners.clear();
            for (int i = 0; i < 4; ++i)
                corners.push_back({uniform(-1.0, 1.0), uniform(-1.0, 1.0), uniform(-1.0, 1.0)});
            if (simplexMeasure(corners[0], corners[1], corners[2], corners[3]) > 0.1)
                break;
        }
        for (auto& x : corners)
            x = place(x);
        return corners;
    }

    // bottom polygon in the plane z = -0.5, top polygon on the lines towards an apex (frustum) or translated
    auto bottom = ellipsePoints(type.isPrism() ? 3 : 4);
    if (!type.isPrism())
        std::swap(bottom[2], bottom[3]);
    for (const auto& q : bottom)
        corners.push_back({q[0], q[1], -0.5});

    const Point<3> apex{uniform(-0.3, 0.3), uniform(-0.3, 0.3), uniform(0.5, 1.5)};
    if (type.isPyramid())
        corners.push_back(apex);
    else
    {
        const bool translated = uniform(0.0, 1.0) < 0.3;
        const double lambda = uniform(0.3, 0.6);
        const Point<3> shift{uniform(-0.3, 0.3), uniform(-0.3, 0.3), uniform(0.5, 1.0)};
        const auto n = corners.size();
        for (std::size_t i = 0; i < n; ++i)
            corners.push_back(translated ? corners[i] + shift : corners[i] + lambda*(apex - corners[i]));
    }

    for (auto& x : corners)
        x = place(x);
    return corners;
}

template<int dim>
Point<dim> randomCenter()
{
    Point<dim> c;
    for (auto& x : c)
        x = uniform(0.35, 0.65);
    return c;
}

////////////////////////////////////////////////////////////////////////////////
// structured grids of cubes or simplices on the unit cube (reference coordinates)
////////////////////////////////////////////////////////////////////////////////

template<int dim>
std::vector<std::pair<Dune::GeometryType, std::vector<Point<dim>>>> unitCubeGrid(int n, bool simplices)
{
    std::vector<std::pair<Dune::GeometryType, std::vector<Point<dim>>>> cells;
    const double h = 1.0/n;
    std::array<int, dim> idx{};
    const int numCells = dim == 2 ? n*n : n*n*n;
    for (int k = 0; k < numCells; ++k)
    {
        int rest = k;
        for (int d = 0; d < dim; ++d) { idx[d] = rest % n; rest /= n; }
        std::vector<Point<dim>> cube;
        for (int c = 0; c < (1 << dim); ++c)
        {
            Point<dim> x;
            for (int d = 0; d < dim; ++d)
                x[d] = h*(idx[d] + ((c >> d) & 1));
            cube.push_back(x);
        }

        if (!simplices)
            cells.push_back({Dune::GeometryTypes::cube(dim), cube});
        else if constexpr (dim == 2)
        {
            cells.push_back({Dune::GeometryTypes::triangle, {cube[0], cube[1], cube[3]}});
            cells.push_back({Dune::GeometryTypes::triangle, {cube[0], cube[3], cube[2]}});
        }
        else
        {
            const int tets[6][4] = {{0,4,5,7}, {0,4,7,6}, {0,2,3,7}, {0,2,7,6}, {0,1,3,7}, {0,1,7,5}};
            for (const auto& t : tets)
                cells.push_back({Dune::GeometryTypes::tetrahedron, {cube[t[0]], cube[t[1]], cube[t[2]], cube[t[3]]}});
        }
    }
    return cells;
}

template<int mydim, int dim, class Map>
std::vector<Geo<mydim, dim>> mapCells(const std::vector<std::pair<Dune::GeometryType, std::vector<Point<mydim>>>>& cells,
                                      const Map& map)
{
    std::vector<Geo<mydim, dim>> geos;
    for (const auto& [type, corners] : cells)
    {
        std::vector<Point<dim>> c;
        for (const auto& p : corners)
            c.push_back(map(p));
        geos.emplace_back(type, c);
    }
    return geos;
}

template<class G>
auto makeTree(const std::vector<G>& geometries)
{
    using Set = Dumux::GeometriesEntitySet<G>;
    return std::make_shared<Dumux::BoundingBoxTree<Set>>(std::make_shared<Set>(geometries));
}

////////////////////////////////////////////////////////////////////////////////
// property: the intersections with all cells of a grid add up to the measure of the geometry
////////////////////////////////////////////////////////////////////////////////

template<class Policy, class G, class Cell>
void checkPartition(const std::string& what, const std::vector<G>& geometries, const std::vector<Cell>& cells)
{
    const auto cellTree = makeTree(cells);
    double maxError = 0.0, maxErrorTrees = 0.0, maxErrorOrder = 0.0;
    for (const auto& g : geometries)
    {
        const double measure = elementMeasure(g);

        // geometry against tree, tree against tree in both orders
        const double sum = sumOfPieces(Dumux::intersectingEntities(g, *cellTree, Policy{}));
        const auto gTree = makeTree(std::vector<G>{g});
        const double sumAB = sumOfPieces(Dumux::intersectingEntities(*gTree, *cellTree, Policy{}));
        const double sumBA = sumOfPieces(Dumux::intersectingEntities(*cellTree, *gTree, Policy{}));
        using std::abs; using std::max;
        maxError = max(maxError, abs(sum - measure)/measure);
        maxErrorTrees = max({maxErrorTrees, abs(sumAB - measure)/measure, abs(sumBA - measure)/measure});

        // both argument orders, cell by cell
        using AB = Dumux::GeometryIntersection<G, Cell, Policy>;
        using BA = Dumux::GeometryIntersection<Cell, G, Policy>;
        for (const auto& cell : cells)
        {
            typename AB::Intersection isAB; typename BA::Intersection isBA;
            const bool foundAB = AB::intersection(g, cell, isAB);
            const bool foundBA = BA::intersection(cell, g, isBA);
            const double mAB = foundAB ? intersectionMeasure<Policy>(isAB) : 0.0;
            const double mBA = foundBA ? intersectionMeasure<Policy>(isBA) : 0.0;
            maxErrorOrder = max(maxErrorOrder, abs(mAB - mBA)/measure);
        }
    }

    report(what + ", geometry and tree", maxError, 1e-10, geometries.size());
    report(what + ", two trees", maxErrorTrees, 1e-10, geometries.size());
    report(what + ", argument order", maxErrorOrder, 1e-10, geometries.size()*cells.size());
}

template<int dim>
void partitionTests(const Transformation<dim>& t)
{
    using namespace Dumux::IntersectionPolicy;
    constexpr int numSamples = 12;
    const std::string where = std::to_string(dim) + "d, " + t.name() + ": ";

    for (const bool simplices : {false, true})
    {
        const std::string grid = simplices ? " / simplex grid" : " / cube grid";
        const auto cells = mapCells<dim, dim>(unitCubeGrid<dim>(4, simplices), t);

        // elements of the world dimension
        std::vector<Dune::GeometryType> types{Dune::GeometryTypes::simplex(dim), Dune::GeometryTypes::cube(dim)};
        if constexpr (dim == 3)
        {
            types.push_back(Dune::GeometryTypes::pyramid);
            types.push_back(Dune::GeometryTypes::prism);
        }
        for (const auto& type : types)
        {
            std::vector<Geo<dim, dim>> geos;
            for (int i = 0; i < numSamples; ++i)
            {
                if constexpr (dim == 2)
                    geos.push_back(makeGeometry<2, 2>(type, randomPolygon<2>(type, randomCenter<2>(), 0.12, randomRotation<2>()), t));
                else
                    geos.push_back(makeGeometry<3, 3>(type, randomPolyhedron(type, randomCenter<3>(), 0.1), t));
            }
            using Policy = std::conditional_t<dim == 2, PolygonPolicy<double, dim>, PolyhedronPolicy<double, dim>>;
            std::ostringstream name; name << where << type << grid;
            checkPartition<Policy>(name.str(), geos, cells);
        }

        // polygons in arbitrary planes in 3d
        if constexpr (dim == 3)
        {
            for (const auto& type : {Dune::GeometryTypes::triangle, Dune::GeometryTypes::quadrilateral})
            {
                std::vector<Geo<2, 3>> geos;
                for (int i = 0; i < numSamples; ++i)
                    geos.push_back(makeGeometry<2, 3>(type, randomPolygon<3>(type, randomCenter<3>(), 0.15, randomRotation<3>()), t));
                std::ostringstream name; name << where << type << " in 3d" << grid;
                checkPartition<PolygonPolicy<double, 3>>(name.str(), geos, cells);
            }
        }

        // segments
        std::vector<Geo<1, dim>> segments;
        for (int i = 0; i < numSamples; ++i)
        {
            const auto c = randomCenter<dim>();
            Point<dim> d;
            for (auto& x : d)
                x = uniform(-0.25, 0.25);
            segments.push_back(makeGeometry<1, dim>(Dune::GeometryTypes::line, {c - d, c + d}, t));
        }
        checkPartition<SegmentPolicy<double, dim>>(where + "segment" + grid, segments, cells);
    }
}

//! polygons and segments in a plane intersected with a grid in the same plane in 3d
void coplanarPartitionTests(const Transformation<3>& t)
{
    using namespace Dumux::IntersectionPolicy;
    constexpr int numSamples = 12;
    const double z = uniform(0.2, 0.8);
    const auto embed = [&](const Point<2>& p){ return t(Point<3>{p[0], p[1], z}); };
    const std::string where = "3d coplanar, " + t.name() + ": ";
    for (const bool simplices : {false, true})
    {
        const std::string grid = simplices ? " / triangle grid" : " / quadrilateral grid";
        const auto cells = mapCells<2, 3>(unitCubeGrid<2>(4, simplices), embed);
        for (const auto& type : {Dune::GeometryTypes::triangle, Dune::GeometryTypes::quadrilateral})
        {
            std::vector<Geo<2, 3>> geos;
            for (int i = 0; i < numSamples; ++i)
                geos.push_back(makeEmbeddedGeometry<2>(type, randomPolygon<2>(type, randomCenter<2>(), 0.12, randomRotation<2>()), embed));
            std::ostringstream name; name << where << type << grid;
            checkPartition<PolygonPolicy<double, 3>>(name.str(), geos, cells);
        }

        std::vector<Geo<1, 3>> segments;
        for (int i = 0; i < numSamples; ++i)
        {
            const auto c = randomCenter<2>();
            const Point<2> d{uniform(-0.25, 0.25), uniform(-0.25, 0.25)};
            segments.push_back(makeEmbeddedGeometry<1>(Dune::GeometryTypes::line, {c - d, c + d}, embed));
        }
        checkPartition<SegmentPolicy<double, 3>>(where + "segment" + grid, segments, cells);
    }
}

////////////////////////////////////////////////////////////////////////////////
// property: point-like intersections of segments with segments and polygons
////////////////////////////////////////////////////////////////////////////////

//! random segment pairs in a plane (z = 0.5 before the transformation in 3d), crossing point from a 2x2 system
template<int dim>
void segmentSegmentPointTests(const Transformation<dim>& t)
{
    using Policy = Dumux::IntersectionPolicy::PointPolicy<double, dim>;
    using Seg = Geo<1, dim>;
    using Algo = Dumux::GeometryIntersection<Seg, Seg, Policy>;
    const auto embed = [&](const Point<2>& p){
        if constexpr (dim == 2) return t(p);
        else return t(Point<3>{p[0], p[1], 0.5});
    };

    std::size_t wrong = 0, samples = 0; double maxError = 0.0;
    for (int k = 0; k < 2000; ++k)
    {
        const Point<2> a{uniform(0, 1), uniform(0, 1)}, b{uniform(0, 1), uniform(0, 1)};
        const Point<2> c{uniform(0, 1), uniform(0, 1)}, d{uniform(0, 1), uniform(0, 1)};
        const Point<2> u = b - a, v = d - c, w = c - a;
        const double det = u[0]*(-v[1]) - u[1]*(-v[0]);
        if (std::abs(det) < 1e-3*u.two_norm()*v.two_norm())
            continue;
        const double s = (w[0]*(-v[1]) - w[1]*(-v[0]))/det;
        const double r = (u[0]*w[1] - u[1]*w[0])/det;
        const double margin = 1e-4;
        const bool inside = s > margin && s < 1 - margin && r > margin && r < 1 - margin;
        const bool outside = s < -margin || s > 1 + margin || r < -margin || r > 1 + margin;
        if (!inside && !outside)
            continue;

        ++samples;
        const Seg seg1(Dune::GeometryTypes::line, std::vector<Point<dim>>{embed(a), embed(b)});
        const Seg seg2(Dune::GeometryTypes::line, std::vector<Point<dim>>{embed(c), embed(d)});
        typename Algo::Intersection p;
        const bool found = Algo::intersection(seg1, seg2, p);
        if (found != inside)
            ++wrong;
        else if (found)
            maxError = std::max(maxError, (p - embed(a + s*u)).two_norm()/(t.scale*u.two_norm()));
    }
    const std::string where = std::to_string(dim) + "d, " + t.name() + ": segment-segment point";
    reportCount(where + " found", wrong, samples);
    report(where + " position", maxError, 1e-10, samples);
}

//! random polygons and segments in 3d, piercing point from the plane of the polygon
void polygonSegmentPointTests(const Transformation<3>& t)
{
    using Policy = Dumux::IntersectionPolicy::PointPolicy<double, 3>;
    std::size_t wrong = 0, samples = 0; double maxError = 0.0;
    for (int k = 0; k < 2000; ++k)
    {
        const auto type = k % 2 ? Dune::GeometryTypes::triangle : Dune::GeometryTypes::quadrilateral;
        const auto corners = randomPolygon<3>(type, randomCenter<3>(), 0.2, randomRotation<3>());
        const Point<3> a{uniform(0.1, 0.9), uniform(0.1, 0.9), uniform(0.1, 0.9)};
        const Point<3> b{uniform(0.1, 0.9), uniform(0.1, 0.9), uniform(0.1, 0.9)};

        // reference: intersection with the plane and signed distances to the edges in the plane
        auto n = Dumux::crossProduct(corners[1] - corners[0], corners[2] - corners[0]);
        n /= n.two_norm();
        const double denom = n*(b - a);
        if (std::abs(denom) < 1e-3*(b - a).two_norm())
            continue;
        const double s = n*(corners[0] - a)/denom;
        const Point<3> p = a + s*(b - a);
        std::vector<Point<3>> cyclic = corners;
        if (type.isQuadrilateral())
            std::swap(cyclic[2], cyclic[3]);
        Point<3> center(0.0);
        for (const auto& x : cyclic)
            center.axpy(1.0/cyclic.size(), x);
        double minDistance = 1e100;
        for (std::size_t i = 0; i < cyclic.size(); ++i)
        {
            const auto e = cyclic[(i + 1) % cyclic.size()] - cyclic[i];
            auto inward = Dumux::crossProduct(n, e);
            inward /= inward.two_norm();
            if (inward*(center - cyclic[i]) < 0)
                inward *= -1.0;
            minDistance = std::min(minDistance, inward*(p - cyclic[i]));
        }
        const double margin = 1e-4;
        const bool inside = s > margin && s < 1 - margin && minDistance > margin*0.2;
        const bool outside = s < -margin || s > 1 + margin || minDistance < -margin*0.2;
        if (!inside && !outside)
            continue;

        ++samples;
        const auto polygon = makeGeometry<2, 3>(type, corners, t);
        const auto segment = makeGeometry<1, 3>(Dune::GeometryTypes::line, {a, b}, t);
        using AB = Dumux::GeometryIntersection<Geo<2, 3>, Geo<1, 3>, Policy>;
        using BA = Dumux::GeometryIntersection<Geo<1, 3>, Geo<2, 3>, Policy>;
        Point<3> isAB, isBA;
        const bool foundAB = AB::intersection(polygon, segment, isAB);
        const bool foundBA = BA::intersection(segment, polygon, isBA);
        if (foundAB != inside || foundBA != inside)
            ++wrong;
        else if (foundAB)
            maxError = std::max({maxError, (isAB - t(p)).two_norm()/(t.scale*(b - a).two_norm()),
                                           (isBA - t(p)).two_norm()/(t.scale*(b - a).two_norm())});
    }
    const std::string where = "3d, " + t.name() + ": polygon-segment point";
    reportCount(where + " found", wrong, samples);
    report(where + " position", maxError, 1e-10, samples);
}

////////////////////////////////////////////////////////////////////////////////
// property: segments touching a polygon (2d) or polyhedron (3d) in a single point
////////////////////////////////////////////////////////////////////////////////

template<int dim>
void contactTests(const Transformation<dim>& t)
{
    using namespace Dumux::IntersectionPolicy;
    using G = Geo<dim, dim>;
    using Seg = Geo<1, dim>;
    using PointAlgo = Dumux::GeometryIntersection<G, Seg, PointPolicy<double, dim>>;
    using SegAlgo = Dumux::GeometryIntersection<G, Seg, SegmentPolicy<double, dim>>;

    std::vector<Dune::GeometryType> types{Dune::GeometryTypes::simplex(dim), Dune::GeometryTypes::cube(dim)};
    if constexpr (dim == 3)
    {
        types.push_back(Dune::GeometryTypes::pyramid);
        types.push_back(Dune::GeometryTypes::prism);
    }

    std::size_t wrongContact = 0, wrongGap = 0, wrongCrossing = 0, wrongSegment = 0, samples = 0;
    double maxError = 0.0;
    for (const auto& type : types)
    {
        for (int k = 0; k < 100; ++k)
        {
            std::vector<Point<dim>> corners;
            if constexpr (dim == 2)
                corners = randomPolygon<2>(type, randomCenter<2>(), 0.2, randomRotation<2>());
            else
                corners = randomPolyhedron(type, randomCenter<3>(), 0.15);
            const G g(type, corners);
            Point<dim> center(0.0);
            for (const auto& x : corners)
                center.axpy(1.0/corners.size(), x);

            // a random point on a random face and the outward normal of that face
            const auto& ref = Dune::referenceElement<double, dim>(type);
            const int face = std::min(int(uniform(0.0, ref.size(1))), ref.size(1) - 1);
            std::vector<Point<dim>> fc;
            for (int j = 0; j < ref.size(face, 1, dim); ++j)
                fc.push_back(corners[ref.subEntity(face, 1, j, dim)]);
            Point<dim> q, normal;
            if constexpr (dim == 2)
            {
                q = fc[0] + uniform(0.1, 0.9)*(fc[1] - fc[0]);
                normal = Point<2>{fc[1][1] - fc[0][1], fc[0][0] - fc[1][0]};
            }
            else
            {
                double l1 = uniform(0.1, 0.8), l2 = uniform(0.1, 0.8);
                if (const double sum = l1 + l2; sum > 0.9) { l1 *= 0.9/sum; l2 *= 0.9/sum; }
                q = fc[0] + l1*(fc[1] - fc[0]) + l2*(fc[2] - fc[0]);
                normal = Dumux::crossProduct(fc[1] - fc[0], fc[2] - fc[0]);
            }
            normal /= normal.two_norm();
            Point<dim> faceCenter(0.0);
            for (const auto& x : fc)
                faceCenter.axpy(1.0/fc.size(), x);
            if (normal*(faceCenter - center) < 0)
                normal *= -1.0;

            // outward direction with a tangential part
            Point<dim> dir = normal;
            for (auto& x : dir)
                x += uniform(-0.5, 0.5);
            if (dir*normal < 0.3*dir.two_norm())
                dir = normal;
            dir *= 0.2/dir.two_norm();

            ++samples;
            const auto seg = [&](const Point<dim>& a, const Point<dim>& b){
                return Seg(Dune::GeometryTypes::line, std::vector<Point<dim>>{t(a), t(b)});
            };
            const G tg = makeGeometry<dim, dim>(type, corners, t);
            typename PointAlgo::Intersection p;
            typename SegAlgo::Intersection s;

            // touching: the contact point for both orientations, nothing with the segment policy
            for (const auto& contact : {seg(q, q + dir), seg(q + dir, q)})
            {
                if (!PointAlgo::intersection(tg, contact, p))
                    ++wrongContact;
                else
                    maxError = std::max(maxError, (p - t(q)).two_norm()/(t.scale*dir.two_norm()));
                if (SegAlgo::intersection(tg, contact, s))
                    ++wrongSegment;
            }
            // a gap: nothing
            if (PointAlgo::intersection(tg, seg(q + 1e-3*dir, q + dir), p))
                ++wrongGap;
            // crossing the boundary: a segment-like intersection, no point
            if (PointAlgo::intersection(tg, seg(q - dir, q + dir), p))
                ++wrongCrossing;
        }
    }
    const std::string where = std::to_string(dim) + "d, " + t.name() + ": segment touching an element";
    reportCount(where + ", contact point found", wrongContact, 2*samples);
    report(where + ", contact point position", maxError, 1e-10, 2*samples);
    reportCount(where + ", no segment-like intersection", wrongSegment, 2*samples);
    reportCount(where + ", gap", wrongGap, samples);
    reportCount(where + ", crossing segment has no point", wrongCrossing, samples);
}

////////////////////////////////////////////////////////////////////////////////
// property: tree queries agree with testing all pairs
////////////////////////////////////////////////////////////////////////////////

template<class Policy, class GA, class GB>
void checkTreeAgainstAllPairs(const std::string& what, const std::vector<GA>& a, const std::vector<GB>& b)
{
    using Algo = Dumux::GeometryIntersection<GA, GB, Policy>;

    // measure (or point) per pair from the tree query
    std::map<std::pair<std::size_t, std::size_t>, double> fromTree, fromPairs;
    std::map<std::pair<std::size_t, std::size_t>, Point<Policy::dimWorld>> pointsTree, pointsPairs;
    for (const auto& is : Dumux::intersectingEntities(*makeTree(a), *makeTree(b), Policy{}))
    {
        fromTree[{is.first(), is.second()}] += simplexMeasure(is.corners());
        if (is.corners().size() == 1)
            pointsTree[{is.first(), is.second()}] = is.corners()[0];
    }
    for (std::size_t i = 0; i < a.size(); ++i)
        for (std::size_t j = 0; j < b.size(); ++j)
            if (typename Algo::Intersection is; Algo::intersection(a[i], b[j], is))
            {
                fromPairs[{i, j}] += intersectionMeasure<Policy>(is);
                if constexpr (Policy::dimIntersection == 0)
                    pointsPairs[{i, j}] = is;
            }

    std::size_t wrong = 0;
    for (const auto& [pair, m] : fromPairs)
        if (!fromTree.count(pair) || std::abs(fromTree.at(pair) - m) > 1e-12*std::max(m, 1e-300))
            ++wrong;
    for (const auto& [pair, m] : fromTree)
        if (!fromPairs.count(pair))
            ++wrong;
    for (const auto& [pair, p] : pointsPairs)
        if (pointsTree.count(pair) && pointsTree.at(pair) != p)
            ++wrong;
    reportCount(what + " (" + std::to_string(fromPairs.size()) + " intersecting pairs)", wrong, a.size()*b.size());
}

template<int dim>
void treeTests(const Transformation<dim>& t)
{
    using namespace Dumux::IntersectionPolicy;
    const std::string where = std::to_string(dim) + "d, " + t.name() + ": tree query vs all pairs, ";
    std::vector<Geo<1, dim>> segA, segB;
    for (int i = 0; i < 40; ++i)
        for (auto* set : {&segA, &segB})
        {
            Point<dim> a, b;
            for (int d = 0; d < dim; ++d) { a[d] = uniform(0, 1); b[d] = uniform(0, 1); }
            if constexpr (dim == 3)
            {
                // segments in a common plane, so that they cross
                a[2] = 0.5; b[2] = 0.5;
            }
            set->push_back(makeGeometry<1, dim>(Dune::GeometryTypes::line, {a, b}, t));
        }
    checkTreeAgainstAllPairs<PointPolicy<double, dim>>(where + "segments, point policy", segA, segB);

    std::vector<Geo<dim, dim>> elemA, elemB;
    for (int i = 0; i < 30; ++i)
        for (auto* set : {&elemA, &elemB})
        {
            const auto type = i % 2 ? Dune::GeometryTypes::simplex(dim) : Dune::GeometryTypes::cube(dim);
            if constexpr (dim == 2)
                set->push_back(makeGeometry<2, 2>(type, randomPolygon<2>(type, randomCenter<2>(), 0.2, randomRotation<2>()), t));
            else
                set->push_back(makeGeometry<3, 3>(type, randomPolyhedron(type, randomCenter<3>(), 0.15), t));
        }
    using Policy = std::conditional_t<dim == 2, PolygonPolicy<double, dim>, PolyhedronPolicy<double, dim>>;
    checkTreeAgainstAllPairs<Policy>(where + "elements, default policy", elemA, elemB);
    checkTreeAgainstAllPairs<SegmentPolicy<double, dim>>(where + "elements and segments, segment policy", elemA, segA);
    checkTreeAgainstAllPairs<PointPolicy<double, dim>>(where + "elements and segments, point policy", elemA, segA);

    if constexpr (dim == 3)
    {
        std::vector<Geo<2, 3>> polygons;
        for (int i = 0; i < 30; ++i)
        {
            const auto type = i % 2 ? Dune::GeometryTypes::triangle : Dune::GeometryTypes::quadrilateral;
            polygons.push_back(makeGeometry<2, 3>(type, randomPolygon<3>(type, randomCenter<3>(), 0.2, randomRotation<3>()), t));
        }
        std::vector<Geo<1, 3>> segments;
        for (int i = 0; i < 40; ++i)
        {
            const Point<3> a{uniform(0, 1), uniform(0, 1), uniform(0, 1)}, b{uniform(0, 1), uniform(0, 1), uniform(0, 1)};
            segments.push_back(makeGeometry<1, 3>(Dune::GeometryTypes::line, {a, b}, t));
        }
        checkTreeAgainstAllPairs<PointPolicy<double, 3>>(where + "polygons and segments, point policy", polygons, segments);
        checkTreeAgainstAllPairs<PolygonPolicy<double, 3>>(where + "polyhedra and polygons, polygon policy", elemA, polygons);
    }
}

////////////////////////////////////////////////////////////////////////////////
// property: cells of a conforming grid touch but do not intersect, moved into each other they overlap
////////////////////////////////////////////////////////////////////////////////

template<int dim>
void neighborTests(const Transformation<dim>& t)
{
    using Policy = std::conditional_t<dim == 2, Dumux::IntersectionPolicy::PolygonPolicy<double, dim>,
                                                Dumux::IntersectionPolicy::PolyhedronPolicy<double, dim>>;
    using G = Geo<dim, dim>;
    using Algo = Dumux::GeometryIntersection<G, G, Policy>;
    const std::string where = std::to_string(dim) + "d, " + t.name() + ": ";
    for (const bool simplices : {false, true})
    {
        const std::string grid = simplices ? "simplex grid" : "cube grid";
        const auto refCells = unitCubeGrid<dim>(3, simplices);
        const auto cells = mapCells<dim, dim>(refCells, t);
        std::size_t wrong = 0, pairs = 0;
        for (std::size_t i = 0; i < cells.size(); ++i)
            for (std::size_t j = i + 1; j < cells.size(); ++j, ++pairs)
                if (typename Algo::Intersection is; Algo::intersection(cells[i], cells[j], is))
                    ++wrong;
        reportCount(where + "distinct cells of a " + grid + " do not intersect", wrong, pairs);

        // the second cell moved by delta towards the first one along the x axis overlaps it by delta times the face
        if (!simplices)
        {
            const double h = 1.0/3.0;
            double maxError = 0.0;
            for (const double delta : {1e-3, 1e-2, 0.5})
            {
                auto moved = refCells[1].second;
                for (auto& x : moved)
                    x[0] -= delta*h;
                const auto second = makeGeometry<dim, dim>(refCells[1].first, moved, t);
                typename Algo::Intersection is;
                const double expected = delta*h*std::pow(h, dim - 1)*std::pow(t.scale, dim);
                const double m = Algo::intersection(cells[0], second, is) ? intersectionMeasure<Policy>(is) : 0.0;
                maxError = std::max(maxError, std::abs(m - expected)/expected);
            }
            report(where + "cubes moved into each other overlap by the moved volume", maxError, 1e-8, 3);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// property: point queries find the cells containing a point
////////////////////////////////////////////////////////////////////////////////

//! signed distances of a point to the faces of a simplex in reference coordinates (positive inside)
template<int dim>
double minBarycentric(const std::vector<Point<dim>>& simplex, const Point<dim>& p)
{
    // solve p = c0 + sum_i l_i (c_i - c0)
    Dune::FieldMatrix<double, dim, dim> A;
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            A[j][i] = simplex[i + 1][j] - simplex[0][j];
    Point<dim> l;
    A.solve(l, p - simplex[0]);
    double l0 = 1.0, minimum = 1e100;
    for (const auto& x : l)
    {
        l0 -= x;
        minimum = std::min(minimum, x);
    }
    return std::min(minimum, l0);
}

template<int dim>
void pointQueryTests(const Transformation<dim>& t)
{
    const std::string where = std::to_string(dim) + "d, " + t.name() + ": ";
    for (const bool simplices : {false, true})
    {
        const std::string grid = simplices ? "simplex grid" : "cube grid";
        const int n = 4;
        const auto refCells = unitCubeGrid<dim>(n, simplices);
        const auto tree = makeTree(mapCells<dim, dim>(refCells, t));

        std::size_t wrong = 0, samples = 0;
        for (int k = 0; k < 500; ++k)
        {
            Point<dim> p;
            for (auto& x : p)
                x = uniform(0.0, 1.0);

            // a point on a face shared by two cells in every second sample
            if (k % 2)
                p[k % dim] = double(1 + int(uniform(0.0, n - 1)))/n;

            std::vector<std::size_t> inside, boundary;
            for (std::size_t c = 0; c < refCells.size(); ++c)
            {
                double distance;
                const auto& corners = refCells[c].second;
                if (simplices)
                    distance = minBarycentric<dim>(corners, p);
                else
                {
                    distance = 1e100;
                    for (int d = 0; d < dim; ++d)
                        distance = std::min({distance, p[d] - corners[0][d], corners.back()[d] - p[d]});
                }
                if (distance > 1e-9)
                    inside.push_back(c);
                else if (distance > -1e-9)
                    boundary.push_back(c);
            }

            ++samples;
            auto found = Dumux::intersectingEntities(t(p), *tree);
            std::sort(found.begin(), found.end());
            auto expected = inside;
            expected.insert(expected.end(), boundary.begin(), boundary.end());
            std::sort(expected.begin(), expected.end());
            if (found != expected)
                ++wrong;
        }
        reportCount(where + "point query in a " + grid, wrong, samples);
    }
}

} // end anonymous namespace

int main()
{
    for (const auto& t : transformations<2>())
    {
        partitionTests<2>(t);
        segmentSegmentPointTests<2>(t);
        contactTests<2>(t);
        treeTests<2>(t);
        neighborTests<2>(t);
        pointQueryTests<2>(t);
    }

    for (const auto& t : transformations<3>())
    {
        partitionTests<3>(t);
        coplanarPartitionTests(t);
        segmentSegmentPointTests<3>(t);
        polygonSegmentPointTests(t);
        contactTests<3>(t);
        treeTests<3>(t);
        neighborTests<3>(t);
        pointQueryTests<3>(t);
    }

    if (numFailures > 0)
    {
        std::cout << numFailures << " checks failed" << std::endl;
        return 1;
    }

    std::cout << "All checks passed" << std::endl;
    return 0;
}
