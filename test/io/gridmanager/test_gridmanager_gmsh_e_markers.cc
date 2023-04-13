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
    // maybe initialize MPI and/or multithreading backend
    Dumux::initialize(argc, argv);

    Dumux::Parameters::init(argc, argv, "test_gridmanager_gmsh_e_markers.input");

    auto name = Dumux::getParam<std::string>("Problem.Name");
    Dumux::GridManagerTests<GRIDTYPE>::testElementMarkers("gmsh", name);

    return 0;
}
