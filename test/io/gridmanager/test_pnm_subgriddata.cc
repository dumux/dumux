// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \brief test for subgrid data wrapper
 */
#include "config.h"

#include <dumux/common/parameters.hh>
#include <dumux/io/grid/porenetwork/gridmanager.hh>
#include <dumux/io/grid/gridmanager_sub.hh>
#include <dumux/io/grid/porenetwork/subgriddata.hh>
#include <dumux/discretization/porenetwork/gridgeometry.hh>

#include "pnmgridutilities.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

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

    PoreNetwork::GridGeometry<double, SubGridManager::Grid::LeafGridView> gridGeometry(subGridManager.grid().leafGridView());
    gridGeometry.update(*gridData);

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
