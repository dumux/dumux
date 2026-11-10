//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test the queries of bounding box trees on a grid, a set of triangles and a polyline
 *
 * A regular polygon split into triangles and a polyline are intersected with a structured grid.
 * The intersection areas add up to the area of every triangle and of the polygon, and the
 * intersection lengths add up to the length of the polyline.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/geometryintersection.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <dumux/geometry/triangulation.hh>
#include <dumux/geometry/volume.hh>

template<class Corners>
double triangleArea(const Corners& corners)
{
    return Dumux::convexPolytopeVolume<2>(Dune::GeometryTypes::triangle,
                                          [&](unsigned int i){ return corners[i]; });
}

void checkClose(double value, double expected, double tolerance, const char* what)
{
    using std::abs;
    if (abs(value - expected) > tolerance)
        DUNE_THROW(Dune::Exception, what << ": " << value << ", expected " << expected);
}

int main(int argc, char* argv[])
{
    using namespace Dumux;
    Dune::MPIHelper::instance(argc, argv);

    constexpr int numCells = 20;
    constexpr double cellArea = 1.0/(numCells*numCells);

    //! [grid tree]
    using Grid = Dune::YaspGrid<2>;
    Grid grid({1.0, 1.0}, {numCells, numCells});

    using GridEntitySet = GridViewGeometricEntitySet<Grid::LeafGridView>;
    const auto gridEntitySet = std::make_shared<GridEntitySet>(grid.leafGridView());
    const BoundingBoxTree<GridEntitySet> gridTree(gridEntitySet);
    //! [grid tree]

    //! [point query]
    const Dune::FieldVector<double, 2> point({0.43, 0.61});
    const std::vector<std::size_t> cells = intersectingEntities(point, gridTree);
    //! [point query]

    const auto expectedCell = intersectingEntityCartesianGrid(point, {0.0, 0.0}, {1.0, 1.0}, {numCells, numCells});
    if (cells.size() != 1 || cells[0] != expectedCell)
        DUNE_THROW(Dune::Exception, "Point query found " << cells.size() << " cells, expected only cell " << expectedCell);

    //! [geometries tree]
    using Point = Dune::FieldVector<double, 2>;
    using Triangle = Dune::MultiLinearGeometry<double, 2, 2>;
    const Point center({0.47, 0.52});
    const double radius = 0.31;
    constexpr int numTriangles = 24;

    std::vector<Triangle> triangles;
    for (int i = 0; i < numTriangles; ++i)
    {
        const double phi0 = 2.0*M_PI*i/numTriangles, phi1 = 2.0*M_PI*(i + 1)/numTriangles;
        const Point a = center + Point({radius*std::cos(phi0), radius*std::sin(phi0)});
        const Point b = center + Point({radius*std::cos(phi1), radius*std::sin(phi1)});
        triangles.emplace_back(Dune::GeometryTypes::triangle, std::vector<Point>{center, a, b});
    }

    using TriangleSet = GeometriesEntitySet<Triangle>;
    const BoundingBoxTree<TriangleSet> triangleTree(std::make_shared<TriangleSet>(triangles));
    //! [geometries tree]

    //! [tree query]
    std::vector<double> coveredArea(gridEntitySet->size(), 0.0);
    std::vector<double> triangleOverlap(numTriangles, 0.0);
    for (const auto& intersection : intersectingEntities(triangleTree, gridTree))
    {
        const double area = triangleArea(intersection.corners());
        triangleOverlap[intersection.first()] += area;
        coveredArea[intersection.second()] += area;
    }
    //! [tree query]

    const double polygonArea = 0.5*numTriangles*radius*radius*std::sin(2.0*M_PI/numTriangles);
    double totalArea = 0.0;
    for (const double area : coveredArea)
    {
        if (area > cellArea*(1.0 + 1e-12))
            DUNE_THROW(Dune::Exception, "Covered area " << area << " exceeds the cell area " << cellArea);
        totalArea += area;
    }
    checkClose(totalArea, polygonArea, 1e-13, "Covered area of all cells");
    for (int i = 0; i < numTriangles; ++i)
        checkClose(triangleOverlap[i], triangles[i].volume(), 1e-14, "Covered area of a triangle");

    //! [policy]
    using Segment = Dune::MultiLinearGeometry<double, 1, 2>;
    const std::vector<Point> polyline{{0.03, 0.11}, {0.37, 0.29}, {0.61, 0.83}, {0.97, 0.71}};
    std::vector<Segment> segments;
    for (std::size_t i = 0; i + 1 < polyline.size(); ++i)
        segments.emplace_back(Dune::GeometryTypes::line, std::vector<Point>{polyline[i], polyline[i+1]});

    using SegmentSet = GeometriesEntitySet<Segment>;
    const BoundingBoxTree<SegmentSet> segmentTree(std::make_shared<SegmentSet>(segments));
    using Policy = IntersectionPolicy::SegmentPolicy<double, 2>;
    const auto pieces = intersectingEntities(segmentTree, gridTree, Policy{});
    //! [policy]

    double polylineLength = 0.0, piecesLength = 0.0;
    for (const auto& segment : segments)
        polylineLength += segment.volume();
    for (const auto& piece : pieces)
        piecesLength += (piece.corners()[1] - piece.corners()[0]).two_norm();
    checkClose(piecesLength, polylineLength, 1e-13, "Length of the polyline pieces");

    //! [geometry intersection]
    using Quadrilateral = Dune::MultiLinearGeometry<double, 2, 2>;
    const Triangle triangle(Dune::GeometryTypes::triangle,
                            std::vector<Point>{{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}});
    const Quadrilateral square(Dune::GeometryTypes::quadrilateral,
                               std::vector<Point>{{0.25, 0.25}, {1.0, 0.25}, {0.25, 1.0}, {1.0, 1.0}});

    using Algorithm = GeometryIntersection<Triangle, Quadrilateral>;
    Algorithm::Intersection polygon;
    double overlap = 0.0;
    if (Algorithm::intersection(triangle, square, polygon))
        for (const auto& simplex : triangulate<2, 2>(polygon))
            overlap += triangleArea(simplex);
    //! [geometry intersection]

    if (polygon.size() != 3)
        DUNE_THROW(Dune::Exception, "Intersection polygon has " << polygon.size() << " corners, expected 3");
    checkClose(overlap, 0.125, 1e-15, "Overlap of triangle and square");

    std::cout << "All tests passed" << std::endl;
    return 0;
}
