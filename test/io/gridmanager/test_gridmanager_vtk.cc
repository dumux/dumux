//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for vtk interface of the grid creator
 */
#include <config.h>
#include <iostream>
#include <algorithm>
#include <tuple>

#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/common/capabilities.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/gridmanager.hh>

#if DUMUX_HAVE_GRIDFORMAT
#include <dumux/io/gridwriter.hh>
#endif

#include "gridmanagertests.hh"

namespace Dumux {

template<class GridView, class GridData>
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>>
getGridData(const GridView& gridView, const GridData& data,
            const std::vector<std::string>& cellDataFieldNames, const std::vector<std::string>& pointDataFieldNames)
{
    std::vector<std::vector<double>> cellData(cellDataFieldNames.size(), std::vector<double>(gridView.size(0)));

    for (const auto& element : elements(gridView))
    {
        const auto eIdx = gridView.indexSet().index(element);
        for (int i = 0; i < cellDataFieldNames.size(); ++i)
            cellData[i][eIdx] = data.getParameter(element, cellDataFieldNames[i]);
    }

    std::vector<std::vector<double>> pointData(pointDataFieldNames.size(), std::vector<double>(gridView.size(GridView::dimension)));

    for (const auto& vertex : vertices(gridView))
    {
        const auto vIdx = gridView.indexSet().index(vertex);
        for (int i = 0; i < pointDataFieldNames.size(); ++i)
            pointData[i][vIdx] = data.getParameter(vertex, pointDataFieldNames[i]);
    }

    return std::make_pair(std::move(cellData), std::move(pointData));
}

template<class Grid>
void testVTKReader(const std::string& gridName)
{
    std::string suffix = gridName;
    std::transform(gridName.begin(), gridName.end(), suffix.begin(), [](unsigned char c){ return std::tolower(c); });

    GridManager<Grid> gridManager;
    gridManager.init(gridName);
    const auto gridData = gridManager.getGridData();
    const auto& gridView = gridManager.grid().leafGridView();

    const auto rank = gridView.comm().rank();
    std::cout << "Rank " << rank << " has " << gridView.size(0) << " elements and " << gridView.size(Grid::dimension) << " vertices." << std::endl;

    auto cellDataFieldNames = getParamFromGroup<std::vector<std::string>>(gridName, "Grid.CellData", std::vector<std::string>{});
    auto pointDataFieldNames = getParamFromGroup<std::vector<std::string>>(gridName, "Grid.PointData", std::vector<std::string>{});
    std::vector<std::vector<double>> cellData, pointData;
    std::tie(cellData, pointData) = getGridData(gridView, *gridData, cellDataFieldNames, pointDataFieldNames);

    if (rank == 0)
    {
        std::cout << "Choosing " << cellDataFieldNames.size() << " cell data fields." << std::endl;
        for (int i = 0; i < cellDataFieldNames.size(); ++i)
            std::cout << "-- field '" << cellDataFieldNames[i] << "' with size " << cellData[i].size() << std::endl;
        std::cout << "Choosing " << pointDataFieldNames.size() << " point data fields." << std::endl;
        for (int i = 0; i < pointDataFieldNames.size(); ++i)
            std::cout << "-- field '" << pointDataFieldNames[i] << "' with size " << pointData[i].size() << std::endl;
    }

    if (gridView.comm().size() > 1)
        suffix += "-parallel";

#if DUMUX_HAVE_GRIDFORMAT
    std::shared_ptr<Dumux::IO::GridWriter<typename Grid::LeafGridView>> gridWriter;
    if constexpr (Dune::Capabilities::isCartesian<Grid>::v)
        gridWriter = std::make_shared<Dumux::IO::GridWriter<typename Grid::LeafGridView>>(IO::Format::vti, gridView);
    else if constexpr (Grid::dimension == 1)
        gridWriter = std::make_shared<Dumux::IO::GridWriter<typename Grid::LeafGridView>>(IO::Format::vtp, gridView);
    else
        gridWriter = std::make_shared<Dumux::IO::GridWriter<typename Grid::LeafGridView>>(IO::Format::vtu, gridView);

    for (int i = 0; i < cellDataFieldNames.size(); ++i)
        gridWriter->setCellField(cellDataFieldNames[i], [&,i=i](const auto& e) { return cellData[i][gridView.indexSet().index(e)]; });
    for (int i = 0; i < pointDataFieldNames.size(); ++i)
        gridWriter->setPointField(pointDataFieldNames[i], [&,i=i](const auto& v) { return pointData[i][gridView.indexSet().index(v)]; });

    gridWriter->write("test-gridmanager-vtk-" + suffix + "-0");
#else
    Dune::VTKWriter<typename Grid::LeafGridView> vtkWriter(gridView);

    for (int i = 0; i < cellDataFieldNames.size(); ++i)
        vtkWriter.addCellData(cellData[i], cellDataFieldNames[i]);
    for (int i = 0; i < pointDataFieldNames.size(); ++i)
        vtkWriter.addVertexData(pointData[i], pointDataFieldNames[i]);

    vtkWriter.write("test-gridmanager-vtk-" + suffix + "-0");
#endif
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    Dumux::Parameters::init(argc, argv, "test_gridmanager_vtk.input");

    using namespace Dumux;

#ifndef GRIDTYPE
    #if UGGRID

        testVTKReader<Dune::UGGrid<2>>("UGGrid");
        return 0;

    #elif FOAMGRID

        testVTKReader<Dune::FoamGrid<1, 3>>("FoamGrid");
        return 0;

    #elif YASPGRID

        testVTKReader<Dune::YaspGrid<2>>("YaspGrid");
        return 0;

    #else
    #warning "Unknown grid type. Skipping test."
        return 77;
    #endif
#else
    auto name = Dumux::getParam<std::string>("Problem.Name");
    Dumux::GridManagerTests<GRIDTYPE>::testElementMarkers<double>("vtu", name + "-element");
    Dumux::GridManagerTests<GRIDTYPE>::testVertexMarkers<double>("vtu", name + "-vertex");
    return 0;
#endif
}
