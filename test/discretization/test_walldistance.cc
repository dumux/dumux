// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Wall distance tests
 *
 */
#include <config.h>

#include <iostream>

#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dune/grid/yaspgrid.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/geometry/distancefield.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/walldistance.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

namespace Dumux {

template<class WD, class GridGeometry>
void testVertex(const GridGeometry& gridGeometry, const std::string& outputName, const std::string& algoName)
{
    const auto& gridView = gridGeometry.gridView();
    Dune::VTKWriter<std::decay_t<decltype(gridView)>> writer(gridView);
    Dune::Timer timer;
    WD wallDistance(gridGeometry, WD::atVertexCenters);
    std::cout << "Update for box with algorithm: " << algoName
                << " took " << timer.elapsed() << " seconds" << std::endl;
    writer.addVertexData(wallDistance.wallDistance(), "wallDistance");

    // We write out wall indices but they are not tested since they
    // are generally not unique (just to make debugging easier)
    std::vector<std::size_t> wallElementIndices(gridGeometry.gridView().size(GridGeometry::GridView::dimension));
    for (int i = 0; i < wallElementIndices.size(); ++i)
        wallElementIndices[i] = wallDistance.wallData()[i].eIdx;
    writer.addVertexData(wallElementIndices, "wallElementIndex");
    writer.write("result_" + algoName + "_" + outputName + "_vertex");
}

template<class WD, class GridGeometry>
void testElement(const GridGeometry& gridGeometry, const std::string& outputName, const std::string& algoName)
{
    const auto& gridView = gridGeometry.gridView();
    Dune::VTKWriter<std::decay_t<decltype(gridView)>> writer(gridView);
    Dune::Timer timer;

    const auto top = gridGeometry.bBoxMax()[GridGeometry::Grid::dimension-1];
    const auto bottom = gridGeometry.bBoxMin()[GridGeometry::Grid::dimension-1];

    // Do not consider scvfs at top or bottom of domain.
    const auto considerScvf = [&](const auto& fvGeometry, const auto& scvf){
        return scvf.ipGlobal()[GridGeometry::Grid::dimension-1] < top - 1e-6 &&
               scvf.ipGlobal()[GridGeometry::Grid::dimension-1] > bottom + 1e-6;
    };

    WD wallDistance(gridGeometry, WD::atElementCenters, considerScvf);
    std::cout << "Update for tpfa with algorithm: " << algoName
            << " took " << timer.elapsed() << " seconds" << std::endl;
    writer.addCellData(wallDistance.wallDistance(), "wallDistance");

    // We write out wall indices but they are not tested since they
    // are generally not unique (just to make debugging easier)
    std::vector<std::size_t> wallElementIndices(gridGeometry.gridView().size(0));
    for (int i = 0; i < wallElementIndices.size(); ++i)
        wallElementIndices[i] = wallDistance.wallData()[i].eIdx;
    writer.addCellData(wallElementIndices, "wallElementIndex");
    writer.write("result_" + algoName + "_" + outputName + "_element");
}

// test extracting wall distance at element or vertex centers
template<class GridGeometry>
void test(const GridGeometry& gridGeometry, const std::string& outputName)
{
    testVertex<WallDistance<GridGeometry>>(gridGeometry, outputName, "aabb_tree");
    testElement<WallDistance<GridGeometry>>(gridGeometry, outputName, "aabb_tree");
}

// run test case for given grid and two discretizations (boundary is described by scvfs)
template<class Grid>
void run(const std::string& testCase)
{
    using GridManager = GridManager<Grid>;
    using GridView = typename Grid::LeafGridView;
    GridManager gridManager;
    gridManager.init(testCase);
    auto suffix = testCase;
    std::transform(suffix.begin(), suffix.end(), suffix.begin(),
        [](unsigned char c){ return std::tolower(c); }
    );
    static constexpr bool enableCache = true;
    {
        std::cout << "Testing " + testCase + " for cctpfa grid geometry" << std::endl;
        using GridGeometry = CCTpfaFVGridGeometry<GridView, enableCache>;
        auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
        test(*gridGeometry, suffix + "_cctpfa" );
    }
    {
        std::cout << "Testing " + testCase + " for box grid geometry" << std::endl;
        using GridGeometry = BoxFVGridGeometry<double, GridView, enableCache>;
        auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
        test(*gridGeometry, suffix + "_box");
    }
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;

    // initialize MPI and multithreading environment
    // finalize is done automatically on exit
    initialize(argc, argv);

    // initialize params
    Parameters::init(argc, argv);

    const auto testCase = getParam<std::string>("TestCase", "2DCube");
    if (testCase == "3DMesh")
#if HAVE_DUNE_ALUGRID
        run<Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>>(testCase);
#else
        std::cerr << "This case needs Dune::ALUGrid to run!" << std::endl;
#endif
    else if (testCase == "3DCube")
        run<Dune::YaspGrid<3>>(testCase);
    else if (testCase == "2DCube")
        run<Dune::YaspGrid<2>>(testCase);
    else
        DUNE_THROW(Dune::Exception, "Unknown test case " << testCase);

    return 0;
}
