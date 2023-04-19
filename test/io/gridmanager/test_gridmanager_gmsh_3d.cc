//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for gmsh interface of the grid creator
 */
#include <config.h>
#include <iostream>

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>

#include "gridmanagertests.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);
    Parameters::init(argc, argv, "test_gridmanager_gmsh_3d.input");
    const auto name = getParam<std::string>("Problem.Name");
    const auto refine = Dumux::getParam<bool>("Grid.Refine", true);
    GridManagerTests<GRIDTYPE>::testBoundaryAndElementMarkers("gmsh", name, refine);

    return 0;
}
