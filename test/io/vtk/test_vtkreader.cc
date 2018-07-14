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
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/alugrid/grid.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/io/vtk/vtureader.hh>

int main(int argc, char** argv) try
{
    Dune::MPIHelper::instance(argc, argv);

    Dumux::VTUReader vtkReader("test-in.vtu");
    using Grid = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>;
    Dumux::VTUReader::Data cellData, pointData;
    auto grid = vtkReader.readGrid<Grid>(cellData, pointData, /*verbose=*/true);

    std::cout << "Successfully read a grid with "
              << grid->leafGridView().size(0) << " elements, "
              << grid->leafGridView().size(Grid::dimension) << " vertices." << std::endl;

    Dune::VTKWriter<Grid::LeafGridView> vtkWriter(grid->leafGridView());
    for (const auto& data : cellData)
        vtkWriter.addCellData(data.second, data.first);
    for (const auto& data : pointData)
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
