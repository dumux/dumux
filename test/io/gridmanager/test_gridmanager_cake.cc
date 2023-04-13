//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for the cake grid manager
 */

#include "config.h"

#include <string>
#include <iostream>

#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/cakegridmanager.hh>

#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#endif

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#ifndef USEUG
#define USEUG false
#endif

// The grid type
#if HAVE_DUNE_UGGRID && USEUG==1
template<int dim>
using Grid = Dune::UGGrid<dim>;
#elif HAVE_DUNE_ALUGRID
template<int dim>
using Grid = Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming>;
#endif

template<int dim>
void testCakeGridManager(const std::string& name)
{
    // using declarations
    using GridManager = typename Dumux::CakeGridManager<Grid<dim>>;
    GridManager gridManager;

    // make the grid
    Dune::Timer timer;
    gridManager.init();

    const auto& gridView = gridManager.grid().leafGridView();
    if (gridView.comm().rank() == 0)
        std::cout << "Constructing " << dim << "-d cake grid with " << gridView.size(0) << " elements took "
                  << timer.elapsed() << " seconds.\n";

    // write the grid to a .vtk file
    Dune::VTKWriter<typename Grid<dim>::LeafGridView> vtkWriter(gridView);
    vtkWriter.write(name);
}

int main(int argc, char** argv)
{
    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    // first read parameters from input file
    Dumux::Parameters::init(argc, argv, "test_gridmanager_cake.input");
    const auto name = Dumux::getParam<std::string>("Grid.Name");

    // test 3-D
    testCakeGridManager<3>("cake-3d-" + name);

    // test 2-D
    testCakeGridManager<2>("cake-2d-" + name);

    return 0;
}
