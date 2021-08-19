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
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/geometry/distancefield.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/box/fvgridgeometry.hh>
#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/discretization/walldistance.hh>
// #include <dumux/freeflow/rans/boundarysearchwalldistance.hh>

// include either UG or ALU grid if one of the is available
#if HAVE_DUNE_UGGRID
#include <dune/grid/uggrid.hh>
#include <dumux/io/grid/cakegridmanager.hh>
#elif HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dumux/io/grid/cakegridmanager.hh>
#endif

template<class Geometry>
using DistanceFieldWithBoundingSpheres = Dumux::DistanceField<Geometry, true>;

template<class Geometry>
using DistanceFieldWithoutBoundingSpheres = Dumux::DistanceField<Geometry, false>;

template<class WD, class GridGeometry>
void testBox(const GridGeometry& gridGeometry, const std::string& paramGroup, const std::string& algoName)
{
    const auto& gridView = gridGeometry.gridView();
    Dune::VTKWriter<std::decay_t<decltype(gridView)>> writer(gridView);
    Dune::Timer timer;
    WD wallDistance(gridGeometry, WD::atVertexCenters);
    std::cout << "Update for box with algorithm: " << algoName
                << " took " << timer.elapsed() << " seconds" << std::endl;
    writer.addVertexData(wallDistance.wallDistance(), "wallDistance");
    writer.write("result_" + algoName + "_" + paramGroup);
}

template<class WD, class GridGeometry>
void testCC(const GridGeometry& gridGeometry, const std::string& paramGroup, const std::string& algoName)
{
    const auto& gridView = gridGeometry.gridView();
    Dune::VTKWriter<std::decay_t<decltype(gridView)>> writer(gridView);
    Dune::Timer timer;

    const auto top = gridGeometry.bBoxMax()[GridGeometry::Grid::dimension-1];
    const auto bottom = gridGeometry.bBoxMin()[GridGeometry::Grid::dimension-1];

    // Do not consider scvfs at top or bottom of domain.
    const auto considerScvf = [&](const auto& scvf)
    {
        return scvf.ipGlobal()[GridGeometry::Grid::dimension-1] < top - 1e-6 &&
               scvf.ipGlobal()[GridGeometry::Grid::dimension-1] > bottom + 1e-6;
    };

    WD wallDistance(gridGeometry, WD::atElementCenters, considerScvf);
    std::cout << "Update for tpfa with algorithm: " << algoName
            << " took " << timer.elapsed() << " seconds" << std::endl;
    writer.addCellData(wallDistance.wallDistance(), "wallDistance");

    // Add corresponding wall element indices. We only do this for ccTpfa because
    // the indices are generally not unique for box.
    // Attention: This is the process-local wall element index!
    std::vector<std::size_t> wallElementIndices(gridGeometry.numDofs());
    for (int dofIdx = 0; dofIdx < wallElementIndices.size(); ++dofIdx)
        wallElementIndices[dofIdx] = wallDistance.wallData()[dofIdx].eIdx;
    writer.addCellData(wallElementIndices, "wallElementIndex");
    writer.write("result_" + algoName + "_" + paramGroup);
}

template<class GridGeometry>
void test(const GridGeometry& gridGeometry, const std::string& paramGroup)
{
    using namespace Dumux;

    using WDB = WallDistance<GridGeometry, DistanceFieldWithBoundingSpheres>;
    using WDBF = WallDistance<GridGeometry, DistanceFieldWithoutBoundingSpheres>;
    using WDAABB = WallDistance<GridGeometry>;

    if constexpr (GridGeometry::discMethod == DiscretizationMethod::box)
    {
        testBox<WDB>(gridGeometry, paramGroup, "bounding_spheres");
        testBox<WDBF>(gridGeometry, paramGroup, "brute_force");
        testBox<WDAABB>(gridGeometry, paramGroup, "aabb_tree");
    }
    else
    {
        testCC<WDB>(gridGeometry, paramGroup, "bounding_spheres");
        testCC<WDBF>(gridGeometry, paramGroup, "brute_force");
        testCC<WDAABB>(gridGeometry, paramGroup, "aabb_tree");
    }
}

int main(int argc, char** argv)
{
    using namespace Dumux;
    using Scalar = double;
    static constexpr bool enableCache = true;

    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    // initialize params
    Dumux::Parameters::init(argc, argv);

    if (getParam<bool>("Test.RunCakeGrid", false))
    {
#if HAVE_DUNE_UGGRID || HAVE_DUNE_ALUGRID
        // check 3D cake grid
        {
#if HAVE_DUNE_UGGRID
            using Grid = Dune::UGGrid<3>;
#else
            using Grid = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>;
#endif
            using GridManager = Dumux::CakeGridManager<Grid>;
            using GridView = typename Grid::LeafGridView;
            GridManager gridManager;
            gridManager.init("3DCake");

            {
                std::cout << "Testing 3D cake cctpfa" << std::endl;
                using GridGeometry = CCTpfaFVGridGeometry<GridView, enableCache>;
                auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
                test(*gridGeometry, "3D_cctpfa_cake");
            }
            {
                std::cout << "Testing 3D cake box" << std::endl;
                using GridGeometry = BoxFVGridGeometry<Scalar, GridView, enableCache>;
                auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
                test(*gridGeometry, "3D_box_cake");
            }
        }
#endif
    return 0;
    }

    // check 2D grid (square)
    {
        using GridManager = Dumux::GridManager<Dune::YaspGrid<2>>;
        using GridView = typename GridManager::Grid::LeafGridView;
        GridManager gridManager;
        gridManager.init("2D");
        {
            std::cout << "Testing 2D box" << std::endl;
            using GridGeometry = BoxFVGridGeometry<Scalar, GridView, enableCache>;
            auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
            test(*gridGeometry, "2D_box");
        }
        {
            std::cout << "Testing 2D cctpfa" << std::endl;
            using GridGeometry = CCTpfaFVGridGeometry<GridView, enableCache>;
            auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
            test(*gridGeometry, "2D_cctpfa");
        }
    }

    // check 3D grid (cube)
    {
        using GridManager = Dumux::GridManager<Dune::YaspGrid<3>>;
        using GridView = typename GridManager::Grid::LeafGridView;
        GridManager gridManager;
        gridManager.init("3D");
        {
            std::cout << "Testing 3D box" << std::endl;
            using GridGeometry = BoxFVGridGeometry<Scalar, GridView, enableCache>;
            auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
            test(*gridGeometry, "3D_box");
        }
        {
            std::cout << "Testing 3D cctpfa" << std::endl;
            using GridGeometry = CCTpfaFVGridGeometry<GridView, enableCache>;
            auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
            test(*gridGeometry, "3D_cctpfa");
        }
    }

    return 0;
} // end main
