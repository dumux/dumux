#include <config.h>
#include <iostream>
#include <algorithm>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/common/timer.hh>
#include <dune/geometry/affinegeometry.hh>

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/dgffoam.hh>
#endif

#include <dumux/common/exceptions.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <test/geometry/writetriangulation.hh>

namespace Dumux {

template<class Grid>
class BBoxTreeTests
{
    using GridView = typename Grid::LeafGridView;
    using Scalar = typename Grid::ctype;
    enum { dimWorld = Grid::dimensionworld };
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using EntitySet = GridViewGeometricEntitySet<GridView, 0>;
    using BoundingBoxTree = Dumux::BoundingBoxTree<EntitySet>;

public:
    int build(const GridView& gv)
    {
        // build a bounding box tree
        tree_ = std::make_shared<BoundingBoxTree>();
        tree_->build(std::make_shared<EntitySet>(gv));
        return 0;
    }

    int intersectPoint(const GlobalPosition& p, std::size_t expectedCollisions)
    {
        std::cout << "Intersect with point ("<< p <<") ";

        Dune::Timer timer;
        const auto entities = intersectingEntities(p, *tree_);

        std::cout << " --> " << entities.size() << " intersection(s) found ("
                  << expectedCollisions << " expected) in " << timer.elapsed() << " seconds.\n";

        if (entities.size() != expectedCollisions)
        {
            std::cerr << "Point intersection failed: Expected "
                      << expectedCollisions << " and got "
                      << entities.size() << "!\n";
            return 1;
        }
        return 0;
    }

    template<class Geometry>
    int intersectGeometry(const Geometry& g, std::size_t expectedCollisions)
    {
        std::cout << "Intersect with geometry type ("<< g.type() <<") ";

        Dune::Timer timer;
        const auto entities = intersectingEntities(g, *tree_);

        std::cout << " --> " << entities.size() << " intersection(s) found ("
                  << expectedCollisions << " expected) in " << timer.elapsed() << " seconds.\n";

        if (entities.size() != expectedCollisions)
        {
            std::cerr << "Geometry intersection failed: Expected "
                      << expectedCollisions << " and got "
                      << entities.size() << "!\n";
            return 1;
        }
        return 0;
    }

    template <class OtherGridView>
    int intersectTree(const Dumux::BoundingBoxTree<GridViewGeometricEntitySet<OtherGridView>>& otherTree,
                      const OtherGridView& otherGridView,
                      std::size_t expectedUniqueIntersections,
                      bool checkTotalIntersections = false,
                      std::size_t expectedIntersections = 0,
                      bool writeVTKPolyData = false,
                      int runIdx = 0)
    {
        Dune::Timer timer;
        const auto intersections = intersectingEntities(*tree_, otherTree);
        std::cout << "Computed " << intersections.size() << " tree intersections in " << timer.elapsed() << std::endl;
        timer.reset();

        if (checkTotalIntersections)
        {
            if (intersections.size() != expectedIntersections)
            {
                std::cerr << "BoundingBoxTree intersection failed: Expected "
                          << expectedIntersections << " (total) and got "
                          << intersections.size() << "!" <<std::endl;
                return 1;
            }
        }

        std::vector<std::vector<std::vector<GlobalPosition>>> map;
        map.resize(otherGridView.size(0));
        std::vector<std::vector<GlobalPosition>> uniqueIntersections;
        for (const auto& is : intersections)
        {
            bool add = true;
            for (const auto& i : map[is.second()])
            {
                add = !is.cornersMatch(i);
                if (!add) break;
            }
            if(add)
            {
                map[is.second()].push_back(is.corners());
                uniqueIntersections.push_back(is.corners());
            }
        }

        if (writeVTKPolyData)
            writeVTKPolyDataTriangle(uniqueIntersections, "unique_" + std::to_string(runIdx));

        std::cout << "Found " << uniqueIntersections.size() << " unique intersections "
                  << "in " << timer.elapsed() << std::endl;

        if (uniqueIntersections.size() != expectedUniqueIntersections)
        {
            std::cerr << "BoundingBoxTree intersection failed: Expected "
                      << expectedUniqueIntersections << " (unique) and got "
                      << uniqueIntersections.size() << "!" <<std::endl;
            return 1;
        }
        return 0;
    }

    template <class GeometryType, class ExpectedIntersectingElements>
    int intersectTree(const Dumux::BoundingBoxTree<GeometriesEntitySet<GeometryType>>& otherTree,
                      const ExpectedIntersectingElements& expectedIntersectingElements)
    {
        Dune::Timer timer;
        const auto intersections = intersectingEntities(otherTree, *tree_);
        std::cout << "Computed " << intersections.size() << " tree intersections in " << timer.elapsed() << std::endl;

        std::unordered_map<std::size_t, std::vector<std::size_t>> sortedResults;

        for (const auto& i : intersections)
        {
            sortedResults[i.first()].push_back(i.second());
            std::sort(sortedResults[i.first()].begin(), sortedResults[i.first()].end());
            sortedResults[i.first()].erase(std::unique(sortedResults[i.first()].begin(), sortedResults[i.first()].end()), sortedResults[i.first()].end());
        }

        for (int i = 0; i < expectedIntersectingElements.size(); ++i)
        {
            if (expectedIntersectingElements.at(i) != sortedResults[i])
            {
                std::cout << "Error for geometry type  " << otherTree.entitySet().entity(i).geometry().type()  << std::endl;
                for (int j = 0; j < expectedIntersectingElements.at(i).size(); ++j)
                    if (expectedIntersectingElements.at(i)[j] != sortedResults[i][j])
                        std::cout << "expected index " << expectedIntersectingElements.at(i)[j] << ", got " << sortedResults[i][j] << std::endl;
                return 1;
            }
        }
        return 0;
    }

private:
    std::shared_ptr<BoundingBoxTree> tree_;
};

} // end namespace Dumux

int main (int argc, char *argv[])
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // Some aliases two type tags for tests using two grids
    constexpr int dimworld = WORLD_DIMENSION;
    using Grid = Dune::YaspGrid<dimworld>;
    using Scalar = Grid::ctype;
    enum { dimWorld = Grid::dimensionworld };
    enum { dim = Grid::dimension };
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    // collect returns to determine exit code
    std::vector<int> returns;
    Dumux::BBoxTreeTests<Grid> test;

    int runIdx = 0;
    for (const auto scaling : {1e10, 1.0, 1e-3, 1e-10})
    {
        std::cout << std::endl
                      << "Testing with scaling = " << scaling << std::endl
                      << "***************************************"
                      << std::endl;

        // create a cube grid
        {
            const GlobalPosition lowerLeft(0.0);
            const GlobalPosition upperRight(1.0*scaling);
            constexpr int numCellsX = 33;
            std::array<unsigned int, dim> elems; elems.fill(numCellsX);
            auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, elems);

            // Dune::VTKWriter<Grid::LeafGridView> vtkWriter(grid->leafGridView());
            // vtkWriter.write("grid_dim" + std::to_string(dimworld), Dune::VTK::ascii);

            // bboxtree tests using one bboxtree and a point
            returns.push_back(test.build(grid->leafGridView()));
            returns.push_back(test.intersectPoint(GlobalPosition(0.0*scaling), 1));
            returns.push_back(test.intersectPoint(GlobalPosition(1e-3*scaling), 1));
            returns.push_back(test.intersectPoint(GlobalPosition(1.0*scaling/numCellsX), 1<<dimworld));
            returns.push_back(test.intersectPoint(GlobalPosition(1.0*scaling), 1));

            // bboxtree tests using one bboxtree and a geometry
            // TODO add more such tests
#if WORLD_DIMENSION == 3
            std::array<GlobalPosition, 3> corners{{{0.0, 0.0, 0.0}, {0.0, 1.0*scaling, 0.0}, {1.0*scaling, 1.0*scaling, 0.0}}};
            Dune::AffineGeometry<double, 2, WORLD_DIMENSION> geometry(Dune::GeometryTypes::simplex(2), corners);
            returns.push_back(test.intersectGeometry(geometry, 2145)); // (33*33/2 - 33/2)*4 + 33
#endif

            // test intersection of grid with 1D geometries (lines)
            if (dimWorld > 1)
            {
                using GeometryType = Dune::MultiLinearGeometry<Scalar, 1, dimWorld>;
                const std::vector<GlobalPosition> cornersDiagonal{lowerLeft, upperRight};
                std::vector<GlobalPosition> cornersVertical{GlobalPosition(0.4*scaling), GlobalPosition(0.4*scaling)};
                cornersVertical[1][dimWorld-1] = 0.6*scaling;
                GeometryType diagonalLine(Dune::GeometryTypes::line, cornersDiagonal);
                GeometryType verticalLine(Dune::GeometryTypes::line, cornersVertical);
                std::vector<GeometryType> geometries{std::move(diagonalLine), std::move(verticalLine)};
                using GeometriesEntitySet = Dumux::GeometriesEntitySet<GeometryType>;
                GeometriesEntitySet entitySet(std::move(geometries));
                Dumux::BoundingBoxTree<GeometriesEntitySet> geometriesTree(std::make_shared<GeometriesEntitySet>(entitySet));
                std::unordered_map<std::size_t, std::vector<std::size_t>> expectedIntersectingElements;

                if (dimWorld == 2)
                {
                    expectedIntersectingElements[0] = {0, 34, 68, 102, 136, 170, 204, 238, 272, 306, 340, 374, 408,
                                                       442, 476, 510, 544, 578, 612, 646, 680, 714, 748, 782, 816,
                                                       850, 884, 918, 952, 986, 1020, 1054, 1088};
                    expectedIntersectingElements[1] = {442, 475, 508, 541, 574, 607, 640};
                    returns.push_back(test.intersectTree(geometriesTree, expectedIntersectingElements));
                }
                else
                {
                    expectedIntersectingElements[0] = {0, 1123, 2246, 3369, 4492, 5615, 6738, 7861, 8984, 10107, 11230,
                                                       12353, 13476, 14599, 15722, 16845, 17968, 19091, 20214, 21337, 22460,
                                                       23583, 24706, 25829, 26952, 28075, 29198, 30321, 31444, 32567,
                                                       33690, 34813, 35936};
                    expectedIntersectingElements[1] = {14599 , 15688, 16777, 17866, 18955, 20044, 21133};
                    returns.push_back(test.intersectTree(geometriesTree, expectedIntersectingElements));
                }
            }

            // test intersection with 2D geometry (triangle)
            if (dimWorld > 1)
            {
                using GeometryType = Dune::MultiLinearGeometry<Scalar, 2, dimWorld>;
                const GlobalPosition corner0 = lowerLeft;
                GlobalPosition corner1 = upperRight; corner1*= 0.1;
                GlobalPosition corner2 = corner1; corner2[dimWorld-1] = lowerLeft[dimWorld-1];
                std::vector<GlobalPosition> cornersTriangle{corner0, corner1, corner2};
                const GeometryType triangle(Dune::GeometryTypes::simplex(2), cornersTriangle);
                using GeometriesEntitySet = Dumux::GeometriesEntitySet<GeometryType>;
                GeometriesEntitySet entitySet({triangle});
                Dumux::BoundingBoxTree<GeometriesEntitySet> geometriesTree(std::make_shared<GeometriesEntitySet>(entitySet));
                std::unordered_map<std::size_t, std::vector<std::size_t>> expectedIntersectingElements;

                if (dimWorld == 2)
                {
                    expectedIntersectingElements[0] = {0, 1, 2, 3, 34, 35, 36, 68, 69, 102};
                    returns.push_back(test.intersectTree(geometriesTree, expectedIntersectingElements));
                }
                else
                {
                    expectedIntersectingElements[0] = {0, 34, 68, 102, 1123, 1157, 1191, 2246, 2280, 3369};
                    returns.push_back(test.intersectTree(geometriesTree, expectedIntersectingElements));
                }
            }

            // test intersection with 2D geometry (quadrilateral)
            if (dimWorld > 1)
            {
                using GeometryType = Dune::AxisAlignedCubeGeometry<Scalar, 2, dimWorld>;
                GlobalPosition lowerLeftCube = upperRight; lowerLeftCube *= 0.4;
                GlobalPosition upperRightCube = lowerLeftCube; upperRightCube[0] += 0.2*scaling; upperRightCube[1] += 0.2*scaling;

                GeometryType cube = [&]()
                {
                    if constexpr (dimWorld == 2)
                        return GeometryType(lowerLeftCube, upperRightCube);
                    else
                    {
                        std::bitset<dimWorld> axes;
                        axes.set(0); axes.set(1);
                        return GeometryType(lowerLeftCube, upperRightCube, axes);
                    }
                }();

                using GeometriesEntitySet = Dumux::GeometriesEntitySet<GeometryType>;
                GeometriesEntitySet entitySet(std::vector<GeometryType>{cube});
                Dumux::BoundingBoxTree<GeometriesEntitySet> geometriesTree(std::make_shared<GeometriesEntitySet>(entitySet));
                std::unordered_map<std::size_t, std::vector<std::size_t>> expectedIntersectingElements;

                if (dimWorld == 2)
                {
                    expectedIntersectingElements[0] = {442, 443, 444, 445, 446, 447, 448, 475, 476, 477, 478, 479, 480, 481,
                                                       508, 509, 510, 511, 512, 513, 514, 541, 542, 543, 544, 545, 546, 547,
                                                       574, 575, 576, 577, 578, 579, 580, 607, 608, 609, 610, 611, 612, 613,
                                                       640, 641, 642, 643, 644, 645, 646};
                    returns.push_back(test.intersectTree(geometriesTree, expectedIntersectingElements));
                }
                else
                {
                    expectedIntersectingElements[0] = {14599, 14600, 14601, 14602, 14603, 14604, 14605, 14632, 14633, 14634,
                                                       14635, 14636, 14637, 14638, 14665, 14666, 14667, 14668, 14669,14670,
                                                       14671, 14698, 14699, 14700, 14701, 14702, 14703, 14704, 14731, 14732,
                                                       14733, 14734, 14735, 14736, 14737, 14764, 14765, 14766, 14767,14768,
                                                       14769, 14770, 14797, 14798, 14799, 14800, 14801, 14802, 14803};
                    returns.push_back(test.intersectTree(geometriesTree, expectedIntersectingElements));
                }
            }
        }

#if HAVE_DUNE_FOAMGRID && WORLD_DIMENSION == 3

        ///////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////
        /// 1D-3D TESTS ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////
        {
            const GlobalPosition lowerLeft(0.0);
            const GlobalPosition upperRight(1.0*scaling);
            constexpr int numCellsX = 10;
            std::array<unsigned int, dim> elems; elems.fill(numCellsX);
            auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, elems);

            // Dune::VTKWriter<Grid::LeafGridView> vtkWriter(grid->leafGridView());
            // vtkWriter.write("grid_dim" + std::to_string(dimworld), Dune::VTK::ascii);

            using NetworkGrid = Dune::FoamGrid<1, dimworld>;
            using NetworkGridView = NetworkGrid::LeafGridView;
            using EntitySet = Dumux::GridViewGeometricEntitySet<NetworkGridView, 0>;
            Dumux::BoundingBoxTree<EntitySet> networkTree;

            {
                std::cout << std::endl
                              << "Intersect with other bounding box tree:" << std::endl
                              << "***************************************"
                              << std::endl;

                // create a network grid from gmsh
                auto networkGrid = std::shared_ptr<NetworkGrid>(Dune::GmshReader<NetworkGrid>::read("network1d.msh", false, false));

                // scaling
                for (const auto& vertex : vertices(networkGrid->leafGridView()))
                {
                    auto newPos = vertex.geometry().corner(0);
                    newPos *= scaling;
                    networkGrid->setPosition(vertex, newPos);
                }

                // Dune::VTKWriter<NetworkGridView> lowDimVtkWriter(networkGrid->leafGridView());
                // lowDimVtkWriter.write("network", Dune::VTK::ascii);

                std::cout << "Constructed 1d network grid with " << networkGrid->leafGridView().size(0) << " elements." << std::endl;

                // build the bulk grid bounding box tree
                returns.push_back(test.build(grid->leafGridView()));

                // build the network grid bounding box tree
                networkTree.build(std::make_shared<EntitySet>(networkGrid->leafGridView()));

                // intersect the two bounding box trees
                returns.push_back(test.intersectTree(networkTree, networkGrid->leafGridView(), 20));
            }
            {
                std::cout << std::endl
                              << "Intersect with other bounding box tree (2):" << std::endl
                              << "*******************************************"
                              << std::endl;

                // construct a line network grid
                const GlobalPosition lowerLeftNW({0.5, 0.5, 0.0});
                const GlobalPosition upperRightNW({0.5, 0.5, 1.0});

                // make the grid (structured interval grid in dimworld space)
                Dune::GridFactory<NetworkGrid> factory;

                constexpr auto geomType = Dune::GeometryTypes::line;

                // create a step vector
                auto step = upperRightNW;
                step -= lowerLeftNW, step /= numCellsX;

                // create the vertices
                auto globalPos = lowerLeftNW;
                for (unsigned int vIdx = 0; vIdx <= numCellsX; vIdx++, globalPos += step)
                    factory.insertVertex(globalPos);

                // create the cells
                for(unsigned int eIdx = 0; eIdx < numCellsX; eIdx++)
                    factory.insertElement(geomType, {eIdx, eIdx+1});

                auto networkGrid = std::shared_ptr<NetworkGrid>(factory.createGrid());

                // scaling
                for (const auto& vertex : vertices(networkGrid->leafGridView()))
                {
                    auto newPos = vertex.geometry().corner(0);
                    newPos *= scaling;
                    networkGrid->setPosition(vertex, newPos);
                }

                std::cout << "Constructed 1d network grid with " << networkGrid->leafGridView().size(0) << " elements." << std::endl;

                // build the bulk grid bounding box tree
                returns.push_back(test.build(grid->leafGridView()));

                // build the network grid bounding box tree
                networkTree.build(std::make_shared<EntitySet>(networkGrid->leafGridView()));

                // intersect the two bounding box trees
                returns.push_back(test.intersectTree(networkTree, networkGrid->leafGridView(), 10, true, 40));
            }
        }

        ///////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////
        /// 2D-3D TESTS ///////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////
        {
            const GlobalPosition lowerLeft(0.0);
            const GlobalPosition upperRight(1.0*scaling);
            constexpr int numCellsX = 10;
            std::array<unsigned int, dim> elems; elems.fill(numCellsX);
            auto grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, elems);

            Dune::VTKWriter<Grid::LeafGridView> vtkWriter(grid->leafGridView());
            vtkWriter.write("grid_dim" + std::to_string(dimworld) + "_" + std::to_string(runIdx), Dune::VTK::ascii);

            using NetworkGrid = Dune::FoamGrid<2, dimworld>;
            using NetworkGridView = NetworkGrid::LeafGridView;
            using EntitySet = Dumux::GridViewGeometricEntitySet<NetworkGridView, 0>;
            Dumux::BoundingBoxTree<EntitySet> networkTree;

            {
                std::cout << std::endl
                              << "Intersect with other bounding box tree:" << std::endl
                              << "***************************************"
                              << std::endl;

                // create a network grid from gmsh
                auto networkGrid = std::shared_ptr<NetworkGrid>(Dune::GmshReader<NetworkGrid>::read("network2d.msh", false, false));

                // scaling
                for (const auto& vertex : vertices(networkGrid->leafGridView()))
                {
                    auto newPos = vertex.geometry().corner(0);
                    newPos *= scaling;
                    networkGrid->setPosition(vertex, newPos);
                }

                Dune::VTKWriter<NetworkGridView> lowDimVtkWriter(networkGrid->leafGridView());
                lowDimVtkWriter.write("network_" + std::to_string(runIdx), Dune::VTK::ascii);

                std::cout << "Constructed 2d network grid with " << networkGrid->leafGridView().size(0) << " elements." << std::endl;

                // build the bulk grid bounding box tree
                returns.push_back(test.build(grid->leafGridView()));

                // build the network grid bounding box tree
                networkTree.build(std::make_shared<EntitySet>(networkGrid->leafGridView()));

                // intersect the two bounding box trees
                returns.push_back(test.intersectTree(networkTree, networkGrid->leafGridView(), 342, true, 342, true, runIdx));
            }
        }
#endif
        ++runIdx;
    }

    std::cout << std::endl;

    // determine the exit code
    if (std::any_of(returns.begin(), returns.end(), [](int i){ return i==1; }))
        return 1;

    return 0;
}
