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
 * \brief Test for gmsh interface of the grid creator
 */
#include <config.h>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dumux/common/parameters.hh>

#include "gridmanagertests.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);
    Parameters::init(argc, argv, "test_gridmanager_gmsh_3d.input");
    const auto name = getParam<std::string>("Problem.Name");
    const auto refine = Dumux::getParam<bool>("Grid.Refine", true);
    GridManagerTests<GRIDTYPE>::testBoundaryAndElementMarkers("gmsh", name, refine);

    return 0;
}
