// graham convex hull test + triangulation
#include <config.h>

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <vector>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/timer.hh>

#include <dumux/common/math.hh>
#include <dumux/geometry/triangulation.hh>
#include <dumux/geometry/grahamconvexhull.hh>

#include "writetriangulation.hh"

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
    std::cout << "Computed triangulation of convex hull with " << convexHullPoints.size()
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

    return 0;
}
