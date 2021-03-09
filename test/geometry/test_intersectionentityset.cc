#include <config.h>

#include <memory>
#include <vector>
#include <numeric>

#include <dune/common/exceptions.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectionentityset.hh>

int main (int argc, char *argv[])
{
    using namespace Dumux;

    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    ///////////////////////////////////////////////////////////////////////
    // Extract the boundary of a given grid and intersect it with the grid
    // Check if all inside entities are found correctly (2d-3d in 3d)
    ///////////////////////////////////////////////////////////////////////
    {
        std::cout << "\nIntersect the boundary of a ball grid (sphere) with the ball grid:\n" << std::endl;

        // construct a grid from file
        using Grid = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
        auto gridFactory = std::make_unique<Dune::GridFactory<Grid>>();
        Dune::GmshReader<Grid>::read(*gridFactory, "ball.msh", /*verbose=*/false, /*boundarySegments=*/false);
        auto grid = std::shared_ptr<Grid>(gridFactory->createGrid());
        const auto leafGridView = grid->leafGridView();

        // and a corresponding entity set
        auto gridEntitySet = std::make_shared<GridViewGeometricEntitySet<Grid::LeafGridView>>(leafGridView);

        // create another entity set with the boundary geometries
        using IGeo = Grid::LeafIntersection::Geometry;
        using Geo = Dune::MultiLinearGeometry<IGeo::ctype, IGeo::mydimension, IGeo::coorddimension>;

        // Convert a geometry wrapper obtained from the Dune grid interface
        // to a real geometry with value type semantics (copyable, storable)
        const auto makeGeometry = [](const IGeo& geo){
            std::vector<Geo::GlobalCoordinate> corners(geo.corners());
            for (unsigned int i = 0; i < geo.corners(); ++i)
                corners[i] = geo.corner(i);
            return Geo(geo.type(), std::move(corners));
        };

        std::vector<Geo> boundaryGeometries;
        std::vector<std::size_t> boundaryIndexToElementIndexMap;
        for (const auto& element : elements(leafGridView))
        {
            for (const auto& intersection : intersections(leafGridView, element))
            {
                if (intersection.boundary())
                {
                    boundaryGeometries.push_back(makeGeometry(intersection.geometry()));
                    boundaryIndexToElementIndexMap.push_back(gridEntitySet->index(element));
                }
            }
        }

        // construct the boundary geometry entity set
        auto geoEntitySet = std::make_shared<GeometriesEntitySet<Geo>>(std::move(boundaryGeometries));

        // intersect the grid with it's boundary geometries
        IntersectionEntitySet<GridViewGeometricEntitySet<Grid::LeafGridView>, GeometriesEntitySet<Geo>> intersectionEntitySet;
        intersectionEntitySet.build(gridEntitySet, geoEntitySet);

        // check if intersections make sense
        for (const auto& intersection : intersections(intersectionEntitySet))
        {
            if (intersection.numDomainNeighbors() != 1)
                DUNE_THROW(Dune::Exception, "Intersections of two conforming geometry sets should "
                            << "have exactly one domain neighbor! Found " << intersection.numDomainNeighbors());
            if (intersection.numTargetNeighbors() != 1)
                DUNE_THROW(Dune::Exception, "Intersections of two conforming geometry sets should "
                            << "have exactly one target neighbor! Found " << intersection.numTargetNeighbors());
            const auto domainNeighbor = intersection.domainEntity(0);
            const auto domainNeighborIndex = gridEntitySet->index(domainNeighbor);
            const auto targetNeighbor = intersection.targetEntity(0);
            const auto targetNeighborIndex = geoEntitySet->index(targetNeighbor);
            if (boundaryIndexToElementIndexMap[targetNeighborIndex] != domainNeighborIndex)
                DUNE_THROW(Dune::Exception, "Found unexpected intersection between element " << domainNeighborIndex
                            << " and boundary segment " << targetNeighborIndex
                            << " (with inside element " << boundaryIndexToElementIndexMap[targetNeighborIndex] << ")");

            const auto volIntersection = intersection.geometry().volume();
            const auto volBoundary = targetNeighbor.geometry().volume();
            if (Dune::FloatCmp::ne(volIntersection, volBoundary))
                DUNE_THROW(Dune::Exception, "Volume of intersection and boundary segment have to be equal for conforming geometry!");
        }

        if (intersectionEntitySet.size() != grid->numBoundarySegments())
            DUNE_THROW(Dune::Exception, "More intersections than boundary segments for conforming boundary intersection set");
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Intersect two polylines with same length and different spacing (1d-1d in 3d)
    ////////////////////////////////////////////////////////////////////////////////
    {
        std::cout << "\nIntersect two polyline with same length and given spacing:\n" << std::endl;

        const auto testPolyLineIntersections = [&](int a, int b)
        {
            std::cout << "-- Test number of segments: " << a << " and " << b << "\n";

            using Geo = Dune::MultiLinearGeometry<double, 1, 3>;
            using Point = Geo::GlobalCoordinate;
            const auto origin = Point({1.0, 1.0, 1.0});

            auto makePolyLine = [](const auto& origin, const int numSegments)
            {
                std::vector<Geo> polyLine;
                const auto step = Point({1.0/numSegments, 1.0/numSegments, 1.0/numSegments});
                for (int i = 0; i < numSegments; ++i)
                {
                    std::vector<Point> corners({Point(origin).axpy(i, step), Point(origin).axpy(i+1, step)});
                    polyLine.emplace_back(Dune::GeometryTypes::line, std::move(corners));
                }
                return polyLine;
            };

            // intersect two polylines
            auto geoSet0 = std::make_shared<GeometriesEntitySet<Geo>>(makePolyLine(origin, a));
            auto geoSet1 = std::make_shared<GeometriesEntitySet<Geo>>(makePolyLine(origin, b));
            IntersectionEntitySet<GeometriesEntitySet<Geo>, GeometriesEntitySet<Geo>> intersectionEntitySet;

            intersectionEntitySet.build(geoSet0, geoSet1);

            // greatest common divisor to compute number of intersections
            const auto numIntersections = b + a - std::gcd(a, b);

            if (intersectionEntitySet.size() != numIntersections)
                DUNE_THROW(Dune::Exception, "Wrong number of line segment intersections."
                            << " Expected " << numIntersections << " got " << intersectionEntitySet.size());
        };

        for (int i = 2; i < 10; ++i)
            for (int j = i; j < 10; ++j)
                testPolyLineIntersections(i, j);
    }

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
