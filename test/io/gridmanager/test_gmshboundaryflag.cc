//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 *
 * \brief Test whether BoundarySegmentIndexFlag works as expected with Box and CCTpfa and ALUGrid
 * \note Alu currently defaults to a boundary flag that works for DGF files only
 */
#include <config.h>
#include <iostream>

#if HAVE_DUNE_ALUGRID

#include <dune/alugrid/grid.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundaryflag.hh>
#include <dumux/io/grid/gridmanager_alu.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include "gmshboundaryflagtest.hh"

namespace Dumux {

// In order to use an alternative BoundaryFlag class, we have to adapt the GridGeometryTraits
template<class GridView>
struct MyBoxGridGeometryTraits : public BoxDefaultGridGeometryTraits<GridView>
{
    struct MyScvfTraits : public BoxDefaultScvfGeometryTraits<GridView>
    { using BoundaryFlag = BoundarySegmentIndexFlag; };

    using SubControlVolumeFace = BoxSubControlVolumeFace<GridView, MyScvfTraits>;
};

template<class GridView>
struct MyCCTpfaGridGeometryTraits : public CCTpfaDefaultGridGeometryTraits<GridView>
{
    struct MyScvfTraits : public CCTpfaDefaultScvfGeometryTraits<GridView>
    { using BoundaryFlag = BoundarySegmentIndexFlag; };

    using SubControlVolumeFace = CCTpfaSubControlVolumeFace<GridView, MyScvfTraits>;
};

} // end namespace Dumux


int main(int argc, char** argv)
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);

    Parameters::init(argc, argv, "test_gmshboundaryflag.input");

    // generate the grid manager and initialize the grid
    using Grid = Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();
    auto gridData = gridManager.getGridData();

    ////////////////////////////////////////////////////////////////////////////////////
    // Test using Box discretization
    ////////////////////////////////////////////////////////////////////////////////////

    // generate FV grid geometry
    using BoxFVGridGeometry = BoxFVGridGeometry<double, typename Grid::LeafGridView,
                                             ENABLE_CACHING,
                                             MyBoxGridGeometryTraits<
                                                 typename Grid::LeafGridView
                                                 >
                                             >;
    auto boxFvGridGeometry = std::make_shared<BoxFVGridGeometry>(leafGridView);

    // run the test
    GmshBoundaryFlagTest<Grid>::testGmshBoundaryFlag<BoxFVGridGeometry>(leafGridView, boxFvGridGeometry, gridData);

    ////////////////////////////////////////////////////////////////////////////////////
    // Test using CCTpfa discretization
    ////////////////////////////////////////////////////////////////////////////////////

    // generate FV grid geometry
    using CCTpfaFVGridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView,
                                             ENABLE_CACHING,
                                             MyCCTpfaGridGeometryTraits<
                                                 typename Grid::LeafGridView
                                                 >
                                             >;
    auto ccTpfaFvGridGeometry = std::make_shared<CCTpfaFVGridGeometry>(leafGridView);

    // run the test
    GmshBoundaryFlagTest<Grid>::testGmshBoundaryFlag<CCTpfaFVGridGeometry>(leafGridView, ccTpfaFvGridGeometry, gridData);

    return 0;
}

#endif /* HAVE_DUNE_ALUGRID */
