/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Test for the vtk reader
 */
#include <config.h>
#include <iostream>
#include <memory>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/io/vtk/vtureader.hh>

int main(int argc, char** argv) try
{
    Dune::MPIHelper::instance(argc, argv);

    Dumux::VTKReader vtkReader("test-in.vtu");
    using Grid = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>;
    Dumux::VTKReader::Data cellData, pointData;
    Dune::GridFactory<Grid> gridFactory;
    auto grid = vtkReader.readGrid(gridFactory, cellData, pointData, /*verbose=*/true);
    const auto& gridView = grid->leafGridView();

    std::cout << "Successfully read a grid with "
              << gridView.size(0) << " elements, "
              << gridView.size(Grid::dimension) << " vertices." << std::endl;

    // to write out we need to reorder as the simple addCellData interface expects MCMG mapper indices
    Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView> elementMapper(gridView, Dune::mcmgElementLayout());
    Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView> vertexMapper(gridView, Dune::mcmgVertexLayout());
    std::vector<std::size_t> elementIndex(gridView.size(0));
    std::vector<std::size_t> vertexIndex(gridView.size(Grid::dimension));
    for (const auto& element : elements(gridView))
    {
        elementIndex[gridFactory.insertionIndex(element)] = elementMapper.index(element);
        for (unsigned int i = 0; i < element.subEntities(Grid::dimension); ++i)
        {
            const auto vertex = element.template subEntity<Grid::dimension>(i);
            vertexIndex[gridFactory.insertionIndex(vertex)] = vertexMapper.index(vertex);
        }
    }

    Dumux::VTKReader::Data reorderedCellData = cellData, reorderedPointData = pointData;
    for (const auto& data : cellData)
    {
        auto& reorderedData = reorderedCellData[data.first];
        for (unsigned int i = 0; i < data.second.size(); ++i)
            reorderedData[elementIndex[i]] = data.second[i];
    }

    for (const auto& data : pointData)
    {
        auto& reorderedData = reorderedPointData[data.first];
        for (unsigned int i = 0; i < data.second.size(); ++i)
            reorderedData[vertexIndex[i]] = data.second[i];
    }

    Dune::VTKWriter<Grid::LeafGridView> vtkWriter(gridView);
    for (const auto& data : reorderedCellData)
        vtkWriter.addCellData(data.second, data.first);
    for (const auto& data : reorderedPointData)
        vtkWriter.addVertexData(data.second, data.first);
    vtkWriter.write("test-out");

    return 0;
}
catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
}
catch (std::exception& e) {
    std::cerr << "stdlib reported error: " << e.what() << std::endl;
    return 2;
}
