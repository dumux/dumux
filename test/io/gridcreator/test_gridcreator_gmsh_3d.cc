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
 * \brief Test for gmsh interface of the grid creator
 */
#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dumux/common/parameters.hh>

#include "gridcreatortests.hh"

int main(int argc, char** argv) try
{
    using namespace Dumux;
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);
    Parameters::init(argc, argv, "test_gridcreator_gmsh_3d.input");
    const auto name = getParam<std::string>("Problem.Name");
    GridCreatorTests<GRIDTYPE>::testBoundaryAndElementMarkers("gmsh", name);

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
