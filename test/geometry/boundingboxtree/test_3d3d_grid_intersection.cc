//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#include <config.h>

#include <iostream>

#include <dune/common/timer.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/float_cmp.hh>
#include <dune/common/exceptions.hh>

#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/type.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#endif

#include <dumux/geometry/boundingboxtree.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>

#include <test/geometry/writetriangulation.hh>

namespace Dumux::Test {

template<class GridView1, class GridView2>
auto computeIntersections(const GridView1& gridView1, const std::string& grid1Name,
                          const GridView2& gridView2, const std::string& grid2Name)
{
    // write to file
    Dune::VTKWriter<GridView1> vtkWriter1(gridView1);
    vtkWriter1.write(grid1Name, Dune::VTK::base64);

    Dune::VTKWriter<GridView2> vtkWriter2(gridView2);
    vtkWriter2.write(grid2Name, Dune::VTK::base64);

    using EntitySet1 = GridViewGeometricEntitySet<GridView1, 0>;
    using EntitySet2 = GridViewGeometricEntitySet<GridView2, 0>;
    auto tree1 = std::make_shared<BoundingBoxTree<EntitySet1>>(std::make_shared<EntitySet1>(gridView1));
    auto tree2 = std::make_shared<BoundingBoxTree<EntitySet2>>(std::make_shared<EntitySet2>(gridView2));

    // compute intersections
    Dune::Timer timer;
    const auto treeIntersections = intersectingEntities(*tree1, *tree2);
    std::cout << "Computed " << treeIntersections.size() << " tree intersections in " << timer.elapsed() << std::endl;
    return treeIntersections;
}

template<class Point, class TreeIntersections>
auto convertIntersections(const TreeIntersections& treeIntersections)
{
    Dune::Timer timer;
    std::vector<std::vector<Point>> intersections;
    intersections.reserve(treeIntersections.size());
    for (const auto& is : treeIntersections)
        intersections.emplace_back(std::vector<Point>(is.corners()));
    std::cout << "Converted to output format in " << timer.elapsed() << " seconds." << std::endl;
    return intersections;
}

template<class Point>
void writeIntersections(const std::vector<std::vector<Point>>& intersections, const std::string& name)
{
    Dune::Timer timer;
    std::cout << "Writing " << intersections.size() << " intersections to file ...";
    writeVTUTetrahedron(intersections, name);
    std::cout << " done ( " << timer.elapsed() << " seconds)." << std::endl;
}

template<class Point>
double computeVolume(const std::vector<std::vector<Point>>& intersections)
{
    using Simplex3Geometry = Dune::AffineGeometry<double, 3, 3>;
    double volume = 0.0;
    for (const auto& tet : intersections)
    {
        const auto tetGeo = Simplex3Geometry(
            Dune::GeometryTypes::simplex(3), std::array<Point, 4>{{ tet[0], tet[1], tet[2], tet[3] }}
        );
        volume += tetGeo.volume();
    }
    return volume;
}

template<class GridView1, class GridView2>
void runCheckVolumeTest(const GridView1& gridView1, const GridView2& gridView2, const std::string& testName, double expectedVolume)
{
    using ctype = typename GridView1::Grid::ctype;
    using Point = Dune::FieldVector<ctype, 3>;

    const auto treeIntersections = computeIntersections(
        gridView1, testName + "_1",
        gridView2, testName + "_2"
    );

    const auto intersections = convertIntersections<Point>(treeIntersections);
    writeIntersections<Point>(intersections, testName + "_intersections");

    const double volume = computeVolume<Point>(intersections);
    std::cout << "Computed total intersection volume of " << volume << " (expected " << expectedVolume << ")" << std::endl;
    // for many intersections the sum is large and the precision is slightly lower
    if (Dune::FloatCmp::ne(expectedVolume, volume, 1e-7*intersections.size()))
        DUNE_THROW(Dune::Exception, "Wrong volume " << volume << ", expected " << expectedVolume
            << " (error: " << (volume-expectedVolume)/expectedVolume*100 << "%)");
}

} // end namespace Dumux::Test

int main (int argc, char *argv[])
{
    // maybe initialize mpi
    Dune::MPIHelper::instance(argc, argv);

    using ctype = double;
    using Point = Dune::FieldVector<ctype, 3>;

    {
        // Some aliases two type tags for tests using two grids
        using Grid = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<ctype, 3>>;

        // make the bulk grid
        const Point lowerLeft(0.0);
        const Point upperRight({8.0, 8.0, 5.0});
        constexpr int numCellsPerDir = 5;
        std::array<unsigned int, 3> elems; elems.fill(numCellsPerDir);
        auto coarseGrid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, elems);
        elems.fill(2*numCellsPerDir);
        auto fineGrid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, elems);

        const ctype expectedVolume = 8.0*8.0*5.0;
        Dumux::Test::runCheckVolumeTest(coarseGrid->leafGridView(), fineGrid->leafGridView(), "structured_coarse_fine", expectedVolume);

        const Point lowerLeft2(0.8);
        const Point upperRight2(1.8);
        elems.fill(3);
        auto smallGrid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft2, upperRight2, elems);
        Dumux::Test::runCheckVolumeTest(coarseGrid->leafGridView(), smallGrid->leafGridView(), "unit_inclusion", 1.0);
    }

#if HAVE_DUNE_ALUGRID
    {
        using ALU = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
        using Yasp = Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<ctype, 3>>;

        Dumux::Parameters::init([](Dune::ParameterTree& params){
            params["Grid.File"] = "../../../test/io/gridmanager/grids/unitcube.msh";
        });
        Dumux::GridManager<ALU> gridManagerAlu;
        gridManagerAlu.init();

        const Point lowerLeft(0.5);
        const Point upperRight(1.5);
        constexpr int numCellsPerDir = 2;
        std::array<unsigned int, 3> elems; elems.fill(numCellsPerDir);
        auto yasp = Dune::StructuredGridFactory<Yasp>::createCubeGrid(lowerLeft, upperRight, elems);

        const ctype expectedVolume = 0.5*0.5*0.5;
        Dumux::Test::runCheckVolumeTest(gridManagerAlu.grid().leafGridView(), yasp->leafGridView(), "unstructured_overlap", expectedVolume);
    }
#endif

    return 0;
}
