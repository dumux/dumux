/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \brief Test for the mmesh grid manager
 */
#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dumux/io/grid/gridmanager_mmesh.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;

    // Initialize MPI, finalize is done automatically on exit.
    Dune::MPIHelper::instance(argc, argv);

    // First set parameters
    Dumux::Parameters::init([](auto& params){
        params["Grid.LowerLeft"] = "-1 -1 -1";
        params["Grid.UpperRight"] = "1 1 1";
    });

    {
        using Grid = Dune::MovingMesh<3>;
        GridManager<Grid> gridManager; gridManager.init();
        Dune::VTKWriter<Grid::LeafGridView> vtkWriter(gridManager.grid().leafGridView());
        vtkWriter.write("mmesh-3d");
    }

    Dumux::Parameters::init([](auto& params){
        params["Grid.File"] = "grids/complex_equi_coarse_tri.msh";
    });

    {
        using Grid = Dune::MovingMesh<2>;
        GridManager<Grid> gridManager; gridManager.init();
        Dune::VTKWriter<Grid::LeafGridView> vtkWriter(gridManager.grid().leafGridView());
        vtkWriter.write("mmesh-2d");
    }

    return 0;
}
