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
 * \brief Wall distance tests
 *
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
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>

#include <dumux/freeflow/rans/wallproperties.hh>
#include <dune/alugrid/grid.hh>


template<class GridGeometry>
void test(const GridGeometry& gridGeometry, const std::string& paramGroup)
{
    using namespace Dumux;

    const auto& gridView = gridGeometry.gridView();
    Dune::VTKWriter<std::decay_t<decltype(gridView)>> writer(gridView);
    Dune::Timer timer;

    BoundarySearchWallProperties<GridGeometry> wallProperties(gridGeometry);
    wallProperties.updateWallDistance();

    if constexpr (GridGeometry::discMethod == DiscretizationMethod::box)
        writer.addVertexData(wallProperties.wallDinstance(), "distance_search");
    else
        writer.addCellData(wallProperties.wallDinstance(), "distance_search");

    std::cout << "Boundary search took " << timer.elapsed() << " seconds" << std::endl;
    timer.reset();

    PoissonWallProperties<GridGeometry> poissonWallProperties(gridGeometry);
    poissonWallProperties.updateWallDistance();

    if constexpr (GridGeometry::discMethod == DiscretizationMethod::box)
        writer.addVertexData(poissonWallProperties.wallDinstance(), "distance_poisson");
    else
        writer.addCellData(poissonWallProperties.wallDinstance(), "distance_poisson");

    std::cout << "Poisson problem took " << timer.elapsed() << " seconds" << std::endl;

    if (wallProperties.wallDinstance().size() != poissonWallProperties.wallDinstance().size())
        DUNE_THROW(Dune::InvalidStateException, "Wrong vector sizes");

    writer.write("result_" + paramGroup);
}

int main(int argc, char** argv)
{
    using namespace Dumux;
    using Scalar = double;
    static constexpr bool enableCache = true;

    // parse command line arguments and input file
    Dumux::Parameters::init(argc, argv);

    {
        using GridManager = Dumux::GridManager<Dune::YaspGrid<2>>;
        using GridView = typename GridManager::Grid::LeafGridView;
        GridManager gridManager;
        gridManager.init("2D");

        {
            std::cout << "Testing 2D box" << std::endl;
            using GridGeometry = BoxFVGridGeometry<Scalar, GridView, enableCache>;
            auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
            gridGeometry->update();
            test(*gridGeometry, "2D_box");
        }
        {
            std::cout << "Testing 2D cctpfa" << std::endl;
            using GridGeometry = CCTpfaFVGridGeometry<GridView, enableCache>;
            auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
            gridGeometry->update();
            test(*gridGeometry, "2D_cctpfa");
        }
    }

    {
        using GridManager = Dumux::GridManager<Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>>;
        using GridView = typename GridManager::Grid::LeafGridView;
        GridManager gridManager;
        gridManager.init("3DTriangle");

        {
            using GridGeometry = BoxFVGridGeometry<Scalar, GridView, enableCache>;
            auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
            gridGeometry->update();
            test(*gridGeometry, "3D_box_tria");
        }
        {
            using GridGeometry = CCTpfaFVGridGeometry<GridView, enableCache>;
            auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
            gridGeometry->update();
            test(*gridGeometry, "3D_cctpfa_tria");
        }

    }

    {
        using GridManager = Dumux::GridManager<Dune::YaspGrid<3>>;
        using GridView = typename GridManager::Grid::LeafGridView;
        GridManager gridManager;
        gridManager.init("3D");

        {
            using GridGeometry = BoxFVGridGeometry<Scalar, GridView, enableCache>;
            auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
            gridGeometry->update();
            test(*gridGeometry, "3D_box");
        }
        {
            using GridGeometry = CCTpfaFVGridGeometry<GridView, enableCache>;
            auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
            gridGeometry->update();
            test(*gridGeometry, "3D_cctpfa");
        }
    }

    return 0;
} // end main
