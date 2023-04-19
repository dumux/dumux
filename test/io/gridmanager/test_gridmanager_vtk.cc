//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
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

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/gridmanager_ug.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

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

    return std::make_pair(cellData, pointData);
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

    auto cellDataFieldNames = getParamFromGroup<std::vector<std::string>>(gridName, "Grid.CellData", std::vector<std::string>{});
    auto pointDataFieldNames = getParamFromGroup<std::vector<std::string>>(gridName, "Grid.PointData", std::vector<std::string>{});
    std::vector<std::vector<double>> cellData, pointData;
    std::tie(cellData, pointData) = getGridData(gridView, *gridData, cellDataFieldNames, pointDataFieldNames);

    Dune::VTKWriter<typename Grid::LeafGridView> vtkWriter(gridView);

    for (int i = 0; i < cellDataFieldNames.size(); ++i)
        vtkWriter.addCellData(cellData[i], cellDataFieldNames[i]);
    for (int i = 0; i < pointDataFieldNames.size(); ++i)
        vtkWriter.addVertexData(pointData[i], pointDataFieldNames[i]);

    vtkWriter.write("test-gridmanager-vtk-" + suffix + "-0");
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    Dumux::Parameters::init(argc, argv, "test_gridmanager_vtk.input");

    using namespace Dumux;

#if UGGRID

    testVTKReader<Dune::UGGrid<2>>("UGGrid");
    return 0;

#elif FOAMGRID

    testVTKReader<Dune::FoamGrid<1, 3>>("FoamGrid");
    return 0;
#else
#warning "Unknown grid type. Skipping test."
    return 77;
#endif
}
