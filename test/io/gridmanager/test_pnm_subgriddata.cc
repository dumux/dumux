// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test for subgrid data wrapper
 */
#include "config.h"

#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/io/grid/porenetwork/subgriddata.hh>
#include <dumux/discretization/porenetwork/gridgeometry.hh>

#include "pnmgridutilities.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);

    // parse command line arguments
    Parameters::init(argc, argv);

    constexpr int dimWorld = 3;
    constexpr int dim = 1;
    using HostGridManager = PoreNetwork::GridManager<dimWorld>;

    HostGridManager hostGridManager;
    hostGridManager.init();
    const auto hostGridData = hostGridManager.getGridData();

    writeToVtk("pnm_host_grid", hostGridManager.grid().leafGridView(), hostGridData);

    auto elementSelector = [&](const auto& element)
    {
        return hostGridData->getParameter(element, "ThroatRegionId") == 0;
    };

    using SubGridManager =  GridManager<Dune::SubGrid<dim, HostGridManager::Grid>>;
    SubGridManager subGridManager;
    subGridManager.init(hostGridManager.grid(), elementSelector);

    using SubGridData = PoreNetwork::SubGridData<HostGridManager::Grid, SubGridManager::Grid>;
    auto gridData = std::make_shared<SubGridData>(subGridManager.grid(), hostGridData);

    PoreNetwork::GridGeometry<double, SubGridManager::Grid::LeafGridView> gridGeometry(subGridManager.grid().leafGridView(), *gridData);

    const auto& params = gridData->parameters(*(elements(subGridManager.grid().leafGridView()).begin()));
    for (const auto& p : params)
        std::cout << p << " ";
    std::cout << std::endl;

    for (const auto& element : elements(subGridManager.grid().leafGridView()))
    {
        if (gridData->getParameter(element, "ThroatRegionId") != 0)
        {
            std::cerr << "Element in wrong subregion" << std::endl;
            return 1;
        }
    }

    writeToVtk("pnm_subgrid_0", subGridManager.grid().leafGridView(), gridData);

    return 0;
}
