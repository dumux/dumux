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
 * \brief Test for the cake grid creator
 */
#include "config.h"
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/grid/cakegridcreator.hh>

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

namespace Dumux
{

namespace Properties
{
NEW_TYPE_TAG(GridCreatorCakeTest);
// Set the grid type
#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(GridCreatorCakeTest, Grid, Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>);
#elif HAVE_UG
SET_TYPE_PROP(GridCreatorCakeTest, Grid, Dune::UGGrid<3>);
#endif
}
}

int main(int argc, char** argv) try
{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    // using declarations
    using TypeTag = TTAG(GridCreatorCakeTest);
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridManager = typename Dumux::CakeGridCreator<Grid>;
    GridManager gridManager;

    // first read parameters from input file
    Dumux::Parameters::init(argc, argv, "test_gridcreator_cake.input");

    // make the grid
    Dune::Timer timer;
    gridManager.init();
    std::cout << "Constructing cake grid with " << gridManager.grid().leafGridView().size(0) << " elements took "
              << timer.elapsed() << " seconds.\n";
    // construct a vtk output writer and attach the boundaryMakers
    Dune::VTKWriter<Grid::LeafGridView> vtkWriter(gridManager.grid().leafGridView());
    vtkWriter.write("cake-00000");

    return 0;
}
catch (Dumux::ParameterException &e) {
    std::cerr << e << ". Abort!\n";
    return 1;
}
catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 3;
}
catch (...) {
    std::cerr << "Unknown exception thrown!\n";
    return 4;
}
