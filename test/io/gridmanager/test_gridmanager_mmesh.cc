//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for the mmesh grid manager
 */
#include <config.h>

#include <dune/grid/io/file/vtk.hh>
#include <dumux/common/initialize.hh>
#include <dumux/io/grid/gridmanager_mmesh.hh>

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);

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
