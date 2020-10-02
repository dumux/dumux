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
 *
 * \brief Test for the cake grid manager
 */

#include<string>

#include "config.h"
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/grid/cakegridmanager.hh>

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#ifndef USEUG
#define USEUG false
#endif

// The grid type
#if HAVE_UG && USEUG==1
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
    std::cout << "Constructing " << dim << "-d cake grid with " << gridManager.grid().leafGridView().size(0) << " elements took "
              << timer.elapsed() << " seconds.\n";
    // construct a vtk output writer and attach the boundaryMakers
    Dune::VTKWriter<typename Grid<dim>::LeafGridView> vtkWriter(gridManager.grid().leafGridView());
    vtkWriter.write(name);
}

int main(int argc, char** argv)
{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    // first read parameters from input file
    Dumux::Parameters::init(argc, argv, "test_gridmanager_cake.input");
    const auto name = Dumux::getParam<std::string>("Grid.Name");

    // test 3-D
    testCakeGridManager<3>("cake-3d-" + name);

    // test 2-D
    testCakeGridManager<2>("cake-2d-" + name);

    return 0;
}
