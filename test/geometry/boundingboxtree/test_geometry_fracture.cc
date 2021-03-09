#include <config.h>

#include <iostream>
#include <algorithm>

#include <dune/common/timer.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/dgffoam.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>
#include <test/geometry/writetriangulation.hh>

int main (int argc, char *argv[])
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    // Some aliases two type tags for tests using two grids
    constexpr int dimworld = 3;
    using Grid = Dune::YaspGrid<dimworld>;
    using ctype = Grid::ctype;
    using Point = Dune::FieldVector<ctype, dimworld>;

    // make the bulk grid
    const Point lowerLeft(0.0);
    const Point upperRight({8.0, 8.0, 5.0});
    constexpr int numCellsX = 10;
    std::array<unsigned int, dimworld> elems; elems.fill(numCellsX);
    auto bulkGrid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, elems);

    // write to file
    using GridView = Grid::LeafGridView;
    Dune::VTKWriter<GridView> vtkWriter(bulkGrid->leafGridView());
    vtkWriter.write("bulk", Dune::VTK::base64);

    // make the bulk bounding box tree
    using BulkEntitySet = Dumux::GridViewGeometricEntitySet<GridView, 0>;
    auto bulkTree = std::make_shared<Dumux::BoundingBoxTree<BulkEntitySet>>();
    bulkTree->build(std::make_shared<BulkEntitySet>(bulkGrid->leafGridView()));

    // make the fracture grid
    using FractureGrid = Dune::FoamGrid<2, dimworld>;
    auto fractureGrid = std::shared_ptr<FractureGrid>(Dune::GmshReader<FractureGrid>::read("fracture.msh", false, false));

    // write to file
    using FractureGridView = FractureGrid::LeafGridView;
    Dune::VTKWriter<FractureGridView> lowDimVtkWriter(fractureGrid->leafGridView());
    lowDimVtkWriter.write("fracture", Dune::VTK::base64);

    // make the fracture bounding box tree
    using FractureEntitySet = Dumux::GridViewGeometricEntitySet<FractureGridView, 0>;
    auto fractureTree = std::make_shared<Dumux::BoundingBoxTree<FractureEntitySet>>();
    fractureTree->build(std::make_shared<FractureEntitySet>(fractureGrid->leafGridView()));

    // compute intersections
    Dune::Timer timer;
    const auto treeIntersections = intersectingEntities(*bulkTree, *fractureTree);
    std::cout << "Computed " << treeIntersections.size() << " tree intersections in " << timer.elapsed() << std::endl;
    timer.reset();

    // convert format
    std::vector<std::vector<Point>> intersections;
    intersections.reserve(treeIntersections.size());
    for (const auto& is : treeIntersections)
        intersections.emplace_back(std::vector<Point>(is.corners()));
    std::cout << "Converted to output format in " << timer.elapsed() << " seconds." << std::endl;
    timer.reset();

    std::cout << "Writing " << intersections.size() << "intersections to file ...";
    Dumux::writeVTKPolyDataTriangle(intersections, "bulk_fracture_intersections");
    std::cout << " done ( " << timer.elapsed() << " seconds)." << std::endl;

    return 0;
}
