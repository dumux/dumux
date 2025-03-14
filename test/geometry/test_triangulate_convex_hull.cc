//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
// graham convex hull test + triangulation
#include <config.h>

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>

#include <dumux/common/math.hh>
#include <dumux/geometry/triangulation.hh>
#include <dumux/geometry/grahamconvexhull.hh>

#include "writetriangulation.hh"

using Point3D = Dune::FieldVector<double, 3>;

double tetrahedronVolume(const Point3D& a, const Point3D& b, const Point3D& c, const Point3D& d)
{
    using std::abs;
    return abs(Dumux::crossProduct(b - a, c - a)*(d - a))/6.0;
}

double hullVolume(const std::vector<Point3D>& points)
{
    double volume = 0.0;
    for (const auto& t : Dumux::triangulate<3, 3>(points))
        volume += tetrahedronVolume(t[0], t[1], t[2], t[3]);
    return volume;
}

void checkHullVolume(const std::vector<Point3D>& points, double expected, const std::string& name)
{
    const auto volume = hullVolume(points);
    using std::abs;
    if (abs(volume - expected) > 1e-8*expected)
        DUNE_THROW(Dune::InvalidStateException,
                   "Wrong hull volume for " << name << ": got " << volume << ", expected " << expected);
}

// hull is the tetrahedron spanned by the first four points, the others are interior
void checkInteriorPoints()
{
    checkHullVolume({{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
                     {0.25, 0.25, 0.1}, {0.1, 0.1, 0.1}},
                    1.0/6.0, "tetrahedron with interior points");

    const std::vector<Point3D> nearVertex{
        {-0.58057106851830498, -0.29775076358325248, 0.47892172882172601},
        {-0.76664480766118182, -0.18561860185460144, -0.51301334052091452},
        {-0.35939177696124769, -0.6827810040821285, 0.45321126331016059},
        {-0.621361647200943, -0.046130559471627697, 0.7836681578432878},
        {-0.58057121064396167, -0.29775129553011892, 0.4789039063440802},
        {-0.76662634241292427, -0.18563034701764078, -0.51293196949162612},
        {-0.58795871204284378, -0.34662007135205208, 0.15458067343399479},
        {-0.5234498642758163, -0.42102332668904252, 0.30207239510020373}};
    checkHullVolume(nearVertex,
                    tetrahedronVolume(nearVertex[0], nearVertex[1], nearVertex[2], nearVertex[3]),
                    "tetrahedron with interior points close to a vertex");
}

// the same polyhedron with its planar face perturbed by the round-off an intersection produces
void checkPerturbedPlanarFace()
{
    const std::vector<Point3D> exact{
        {0.0, 0.0, 0.0},
        {0.10000000000000001, 0.20000000000000001, 0.90000000000000002},
        {0.1309712721100709, 0.47874144899063814, 0.71417236733957468},
        {0.1561289667403154, 0.85870931707173481, 0.2341934501104731},
        {0.37130348780858424, 0.94426955653765066, 0.28442695565376508},
        {1.3, 0.10000000000000001, 0.20000000000000001}};

    const std::vector<Point3D> perturbed{
        {0.0, 0.0, 0.0},
        {0.10000000000000001, 0.20000000000000001, 0.90000000000000002},
        {0.13097127211007112, 0.47874144899044035, 0.71417236733980216},
        {0.1561289667403041, 0.85870931707170495, 0.23419345011046996},
        {0.37130348780854766, 0.94426955653761602, 0.28442695565373782},
        {1.3, 0.10000000000000001, 0.20000000000000001}};

    checkHullVolume(perturbed, hullVolume(exact), "polyhedron with a perturbed planar face");
}

// random tetrahedra cut by a random plane, clipped volume known in closed form
void checkClippedTetrahedra()
{
    std::mt19937 gen(20260723);
    std::uniform_real_distribution<double> coord(-1.0, 1.0), frac(0.0, 1.0);

    for (int test = 0; test < 500; ++test)
    {
        std::array<Point3D, 4> corners;
        for (auto& c : corners)
            c = {coord(gen), coord(gen), coord(gen)};
        const auto volume = tetrahedronVolume(corners[0], corners[1], corners[2], corners[3]);
        if (volume < 1e-2)
            continue;

        Point3D normal{coord(gen), coord(gen), coord(gen)};
        if (normal.two_norm() < 0.2)
            continue;
        normal /= normal.two_norm();

        auto lowest = 1e100, highest = -1e100;
        for (const auto& c : corners)
        {
            using std::min; using std::max;
            lowest = min(lowest, c*normal);
            highest = max(highest, c*normal);
        }
        const auto offset = lowest + (highest - lowest)*(0.15 + 0.7*frac(gen));

        std::array<double, 4> distance;
        std::vector<int> inside, outside;
        for (int i = 0; i < 4; ++i)
        {
            distance[i] = corners[i]*normal - offset;
            (distance[i] <= 0.0 ? inside : outside).push_back(i);
        }
        if (inside.empty() || outside.empty())
            continue;

        const auto cutPoint = [&](int i, int j)
        {
            auto p = corners[i];
            p.axpy(distance[i]/(distance[i] - distance[j]), corners[j] - corners[i]);
            return p;
        };
        const auto cornerFraction = [&](int v)
        {
            auto f = 1.0;
            for (int i = 0; i < 4; ++i)
                if (i != v)
                    f *= distance[v]/(distance[v] - distance[i]);
            return f;
        };

        double expected;
        if (inside.size() == 1)
            expected = volume*cornerFraction(inside[0]);
        else if (outside.size() == 1)
            expected = volume*(1.0 - cornerFraction(outside[0]));
        else
        {
            const auto a = inside[0], b = inside[1], p = outside[0], q = outside[1];
            const auto ap = cutPoint(a, p), aq = cutPoint(a, q);
            const auto bp = cutPoint(b, p), bq = cutPoint(b, q);
            expected = tetrahedronVolume(corners[a], corners[b], ap, aq)
                     + tetrahedronVolume(corners[b], ap, aq, bq)
                     + tetrahedronVolume(corners[b], ap, bq, bp);
        }
        if (expected < 1e-3)
            continue;

        std::vector<Point3D> clipped;
        for (int i = 0; i < 4; ++i)
            if (distance[i] <= 0.0)
                clipped.push_back(corners[i]);
        for (int i = 0; i < 4; ++i)
            for (int j = i+1; j < 4; ++j)
                if ((distance[i] > 0.0) != (distance[j] > 0.0))
                    clipped.push_back(cutPoint(i, j));

        const auto id = std::to_string(test);
        checkHullVolume(clipped, expected, "clipped tetrahedron " + id);

        std::shuffle(clipped.begin(), clipped.end(), gen);
        checkHullVolume(clipped, expected, "shuffled clipped tetrahedron " + id);

        Point3D centre(0.0);
        for (const auto& p : clipped)
            centre += p;
        centre /= clipped.size();

        auto withInterior = clipped;
        auto interior = clipped[0];
        interior.axpy(1e-4, centre - clipped[0]);
        withInterior.push_back(interior);
        checkHullVolume(withInterior, expected, "clipped tetrahedron with interior point " + id);

        auto jittered = clipped;
        for (auto& p : jittered)
            for (int i = 0; i < 3; ++i)
                p[i] += 1e-13*coord(gen);
        checkHullVolume(jittered, expected, "jittered clipped tetrahedron " + id);
    }
}

template<class Scalar>
class UniformDistributedRandomNumber
{
public:
    UniformDistributedRandomNumber(Scalar min, Scalar max)
    : gen_(r_()), dist_(min, max) {}

    Scalar operator()()
    { return dist_(gen_); }

private:
    std::random_device r_;
    std::default_random_engine gen_;
    std::uniform_real_distribution<Scalar> dist_;
};

Dune::FieldVector<double, 3>
randomPointOnPlane(const Dune::FieldVector<double, 3>& origin,
                   const Dune::FieldVector<double, 3>& normal,
                   double min = -1.0, double max = 1.0)
{
    using Point = Dune::FieldVector<double, 3>;

    UniformDistributedRandomNumber<double> dice(min, max);
    Point p{dice(), dice(), dice()};
    p = Dumux::crossProduct(normal, p);
    p += origin;
    return p;
}

void writeCSV(const std::vector<Dune::FieldVector<double, 3>>& p,
              const std::string& filename)
{
    std::ofstream fout(filename + ".csv");
    fout << "X, Y, Z\n";
    for (const auto& point : p)
        fout << point[0] << ", " << point[1] << ", " << point[2]  << '\n';
    fout << std::endl;
}

void writeVTKPolyData(const std::vector<Dune::FieldVector<double, 3>>& p,
                      const std::string& filename)
{
    std::ofstream fout(filename + ".vtp");
    fout << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
         << "  <PolyData>\n"
         << "    <Piece NumberOfPoints=\"" << p.size() << "\" NumberOfLines=\"" << p.size() << "\">\n"
         << "      <Points>\n"
         << "        <DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const auto& point : p)
        fout << point << " ";
    fout << '\n';

    fout << "        </DataArray>\n"
         << "      </Points>\n"
         << "      <Lines>\n"
         << "        <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";

    for (std::size_t i = 0; i < p.size()-1; ++i)
        fout << i << " " << i+1 << " ";
    fout << p.size()-1 << " 0\n";

    fout << "        </DataArray>\n";
    fout << "        <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";

    for (std::size_t i = 1; i <= p.size(); ++i)
        fout << i*2 << " ";
    fout << '\n';

    fout << "        </DataArray>\n"
         << "      </Lines>\n"
         << "    </Piece>\n"
         << "</PolyData>\n"
         << "</VTKFile>\n";
}

int main(int argc, char* argv[])
{
    using namespace Dumux;

    using Point = Dune::FieldVector<double, 3>;
    // create 100 random points
    std::vector<Point> points(100);
    // normal vector of the plane
    const Point normal({0.0, 1.0, 1.0});
    const Point origin({3.0, 0.0, 0.0});
    for (auto&& p : points)
        p = randomPointOnPlane(origin, normal, -10.0, 10.0);

    writeCSV(points, "points");

    Dune::Timer timer;
    auto convexHullPoints = grahamConvexHull<2>(points);
    std::cout << "Computed convex hull of " << points.size() << " points in "
              << timer.elapsed() << " seconds." << std::endl;

    writeVTKPolyData(convexHullPoints, "convexhull");

    Dune::Timer timer2;
    auto triangles = triangulate<2, 3>(convexHullPoints);
    std::cout << "Computed 2D (in 3D) triangulation of convex hull with " << convexHullPoints.size()
              << " points in " << timer2.elapsed() << " seconds." << std::endl;

    writeVTKPolyDataTriangle(triangles, "triangulation");

    // some specific tests for corner cases with colinear points
    points = {{0.0,0.0,0.0}, {0.0,0.0,1.0}, {0.0,0.0,2.0}, {0.0,0.0,3.0}};
    auto hull = grahamConvexHull<2>(points);
    if (!(hull.empty())) DUNE_THROW(Dune::InvalidStateException, "False positive for finding a convex hull!");

    points = {{0.0,0.0,0.0}, {0.0,0.0,1.0}, {0.0,0.0,2.0}, {0.0,0.0,3.0}, {0.0,0.0,4.0}, {2.0,3.0,3.0}};
    hull = grahamConvexHull<2>(points);
    if (hull.empty()) DUNE_THROW(Dune::InvalidStateException, "Didn't find convex hull!");

    points = {{2.0,3.0,3.0}, {0.0,0.0,4.0}, {0.0,0.0,0.0}, {0.0,0.0,1.0}, {0.0,0.0,2.0}, {0.0,0.0,3.0}};
    hull = grahamConvexHull<2>(points);
    if (hull.empty()) DUNE_THROW(Dune::InvalidStateException, "Didn't find convex hull!");

    // 3d triangulation tests
    const auto tetPoints = std::vector<Point>{{0.0,0.0,0.0}, {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0}};
    const auto hull3D1 = triangulate<3, 3>(tetPoints);
    writeVTUTetrahedron(hull3D1, "hulltriangulation3d_simplex");
    std::cout << "Created triangulation with " << hull3D1.size() << " tetrahedra" << std::endl;
    if (hull3D1.size() != 1)
        DUNE_THROW(Dune::InvalidStateException, "Incorrect hull for simplex!");

    // in the general case with more than exactly four points, we find triangulations with the mid point rule
    // therefore, we use more than one tetrahedron even if the convex hull is a tetrahedron
    const auto tetPointsInclusion = std::vector<Point>{{0.0,0.0,0.0}, {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0}, {0.25,0.25,0.1}};
    const auto hull3D2 = triangulate<3, 3>(tetPointsInclusion);
    writeVTUTetrahedron(hull3D2, "hulltriangulation3d_simplexincl");
    std::cout << "Created triangulation with " << hull3D2.size() << " tetrahedra" << std::endl;
    if (hull3D2.size() != 4)
        DUNE_THROW(Dune::InvalidStateException, "Incorrect hull for simplex with inclusion!");

    const auto tetPointsCoplanar = std::vector<Point>{{0.0,0.0,0.0}, {1.0,0.0,0.0}, {0.0,1.0,0.0}, {0.0,0.0,1.0}, {0.1,0.1,0.0}};
    const auto hull3D3 = triangulate<3, 3>(tetPointsCoplanar);
    writeVTUTetrahedron(hull3D3, "hulltriangulation3d_coplanar");
    std::cout << "Created triangulation with " << hull3D3.size() << " tetrahedra" << std::endl;
    if (hull3D3.size() != 4)
        DUNE_THROW(Dune::InvalidStateException, "Incorrect hull for simplex with coplanar face inclusion!");

    const auto complexPyramidPoints = std::vector<Point>{
        {0.0,0.0,0.0}, {0.5,0.87,0.0}, {1.5,0.87,0.0},
        {2.0,0.0,0.0}, {1.5,-0.87,0.0}, {0.5,-0.87,0.0},
        {1.0,0.0,1.0}
    };

    const auto hullPyramid = triangulate<3, 3>(complexPyramidPoints);
    writeVTUTetrahedron(hullPyramid, "hulltriangulation3d_pyramid");
    std::cout << "Created triangulation with " << hullPyramid.size() << " tetrahedra" << std::endl;
    if (hullPyramid.size() != 6*2)
        DUNE_THROW(Dune::InvalidStateException, "Incorrect hull for pyramid with hexagonal base!");

    points.clear();
    UniformDistributedRandomNumber<double> dice(-1.0, 1.0);
    for (int i = 0; i < 10; ++i)
        points.emplace_back(Point{ dice(), dice(), dice() });

    Dune::Timer timer3;
    const auto hull3D4Rnd = triangulate<3, 3>(points);
    std::cout << "Created triangulation with " << hull3D4Rnd.size() << " tetrahedra" << std::endl;
    std::cout << "Computed 3D (in 3D) convex hull based triangulation of point cloud with " << points.size()
              << " points in " << timer3.elapsed() << " seconds." << std::endl;
    writeVTUTetrahedron(hull3D4Rnd, "hulltriangulation3d_rnd");
    if (hull3D4Rnd.size() < 4)
        DUNE_THROW(Dune::InvalidStateException, "Incorrect hull for random points!");

    checkInteriorPoints();
    checkPerturbedPlanarFace();
    checkClippedTetrahedra();

    return 0;
}
