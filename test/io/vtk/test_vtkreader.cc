//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for the vtk reader
 */
#include <config.h>
#include <iostream>
#include <memory>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/grid/yaspgrid.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#endif

#include <dumux/common/initialize.hh>
#include <dumux/common/exceptions.hh>
#include <dumux/io/vtk/vtkreader.hh>

int main(int argc, char** argv)
{
    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    if (argc != 3)
        DUNE_THROW(Dune::IOError, "Needs two arguments, the vtk file name and an output file base name");

    auto vtkReader = std::make_shared<Dumux::VTKReader>(std::string(argv[1]));
    using Grid = GRIDTYPE;
    Dumux::VTKReader::Data cellData, pointData;
    Dune::GridFactory<Grid> gridFactory;
    auto grid = vtkReader->readGrid(gridFactory, cellData, pointData, /*verbose=*/true);
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
    for (const auto& [name, data] : cellData)
    {
        // sanity check
        if (!vtkReader->hasData(name, Dumux::VTKReader::DataType::cellData))
            DUNE_THROW(Dune::Exception, "Array " << name << " exists but hasData returns false!");

        auto& reorderedData = reorderedCellData[name];
        for (unsigned int i = 0; i < data.size(); ++i)
            reorderedData[elementIndex[i]] = data[i];
    }

    for (const auto& [name, data] : pointData)
    {
        // sanity check
        if (!vtkReader->hasData(name, Dumux::VTKReader::DataType::pointData))
            DUNE_THROW(Dune::Exception, "Array " << name << " exists but hasData returns false!");

        auto& reorderedData = reorderedPointData[name];
        for (unsigned int i = 0; i < data.size(); ++i)
            reorderedData[vertexIndex[i]] = data[i];
    }

    Dune::VTKWriter<Grid::LeafGridView> vtkWriter(gridView);
    for (const auto& [name, data] : reorderedCellData)
        vtkWriter.addCellData(data, name);
    for (const auto& [name, data] : reorderedPointData)
        vtkWriter.addVertexData(data, name);
    vtkWriter.write(std::string(argv[2]));

    return 0;
}
