#include <config.h>
#include <iostream>
#include <algorithm>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/common/timer.hh>

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
#include <dune/foamgrid/dgffoam.hh>
#else
#include <dune/foamgrid/dgffoam.cc>
#endif
#endif

#include <dumux/common/exceptions.hh>
#include <dumux/common/geometry/boundingboxtree.hh>
#include <dumux/common/geometry/geometricentityset.hh>
#include <dumux/common/geometry/intersectingentities.hh>
#include <test/common/geometry/writetriangulation.hh>

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

    template <class OtherEntitySet, class OtherGridView>
    int intersectTree(const Dumux::BoundingBoxTree<OtherEntitySet>& otherTree,
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

private:
    std::shared_ptr<BoundingBoxTree> tree_;
};

} // end namespace Dumux

int main (int argc, char *argv[]) try
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

            // bboxtree tests using one bboxtree
            returns.push_back(test.build(grid->leafGridView()));
            returns.push_back(test.intersectPoint(GlobalPosition(0.0*scaling), 1));
            returns.push_back(test.intersectPoint(GlobalPosition(1e-3*scaling), 1));
            returns.push_back(test.intersectPoint(GlobalPosition(1.0*scaling/numCellsX), 1<<dimworld));
            returns.push_back(test.intersectPoint(GlobalPosition(1.0*scaling), 1));
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

        #if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
                constexpr auto geomType = Dune::GeometryTypes::line;
        #else
                auto geomType = Dune::GeometryType(1);
        #endif

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
// //////////////////////////////////
//   Error handler
// /////////////////////////////////
catch (Dumux::ParameterException &e) {
    std::cerr << e << ". Abort!\n";
    return 1;
}
catch (const Dune::Exception& e) {
    std::cout << e << std::endl;
    return 1;
}
