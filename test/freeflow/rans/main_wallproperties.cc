// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
 * \ingroup RANSTests
 * \brief Pipe flow test for the staggered grid RANS model,
 *
 * This test simulates is based on pipe flow experiments by
 * John Laufers experiments in 1954 \cite Laufer1954a.
 */
#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>

#include <dumux/freeflow/rans/wallproperties.hh>

// #include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    GridManager<Dune::YaspGrid<2>> gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    static constexpr bool enableCache = true;
    using GridView = std::decay_t<decltype(leafGridView)>;
    using Scalar = double;
    using GridGeometry = BoxFVGridGeometry<Scalar, GridView, enableCache>;

    // create the finite volume grid geometry
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);
    gridGeometry->update();

    // stop time for the entire computation
    Dune::Timer timer;

    BoundarySearchWallProperties<GridGeometry> wallProperties(*gridGeometry);

    wallProperties.updateWallDistance();

    const auto& d = wallProperties.wallDinstance();

    Dune::VTKWriter<GridView> writer(leafGridView);
    writer.addVertexData(d, "distance");

    writer.write("result");

    std::cout << "boundary search took " << timer.elapsed() << " seconds" << std::endl;

    timer.reset();

    PoissonWallProperties<GridGeometry> poissonWallProperties(*gridGeometry);

    poissonWallProperties.solve();

    std::cout << "Poisson problem took " << timer.elapsed() << " seconds" << std::endl;


    return 0;
} // end main
