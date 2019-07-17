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
 * \brief Test whether GmshBoundaryFlag works as expected with Box and CCTpfa
 */
#include <config.h>
#include <iostream>

#if HAVE_DUNE_ALUGRID

#include <dune/alugrid/grid.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/gridmanager.hh>

#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include "gmshboundaryflagtest.hh"

namespace Dumux {


// In order to use an alternative BoundaryFlag class, we have to adapt the GridGeometryTraits
template<class GridView, class MyBoundaryFlag>
struct MyBoxGridGeometryTraits : public BoxDefaultGridGeometryTraits<GridView>
{
    struct MyScvfTraits : public BoxDefaultScvfGeometryTraits<GridView>
    { using BoundaryFlag = MyBoundaryFlag; };

    using SubControlVolumeFace = BoxSubControlVolumeFace<GridView, MyScvfTraits>;
};

template<class GridView, class MyBoundaryFlag>
struct MyCCTpfaGridGeometryTraits : public CCTpfaDefaultGridGeometryTraits<GridView>
{
    struct MyScvfTraits : public CCTpfaDefaultScvfGeometryTraits<GridView>
    { using BoundaryFlag = MyBoundaryFlag; };

    using SubControlVolumeFace = CCTpfaSubControlVolumeFace<GridView, MyScvfTraits>;
};

} // end namespace Dumux


int main(int argc, char** argv) try
{
    using namespace Dumux;

    Dune::MPIHelper::instance(argc, argv);

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
                                                 typename Grid::LeafGridView, GmshBoundaryFlag<Grid>
                                                 >
                                             >;
    auto boxFvGridGeometry = std::make_shared<BoxFVGridGeometry>(leafGridView);
    boxFvGridGeometry->update();

    // run the test
    GmshBoundaryFlagTest<Grid>::testGmshBoundaryFlag<BoxFVGridGeometry>(leafGridView, boxFvGridGeometry, gridData);

    ////////////////////////////////////////////////////////////////////////////////////
    // Test using CCTpfa discretization
    ////////////////////////////////////////////////////////////////////////////////////

    // generate FV grid geometry
    using CCTpfaFVGridGeometry = CCTpfaFVGridGeometry<typename Grid::LeafGridView,
                                             ENABLE_CACHING,
                                             MyCCTpfaGridGeometryTraits<
                                                 typename Grid::LeafGridView, GmshBoundaryFlag<Grid>
                                                 >
                                             >;
    auto ccTpfaFvGridGeometry = std::make_shared<CCTpfaFVGridGeometry>(leafGridView);
    ccTpfaFvGridGeometry->update();

    // run the test
    GmshBoundaryFlagTest<Grid>::testGmshBoundaryFlag<CCTpfaFVGridGeometry>(leafGridView, ccTpfaFvGridGeometry, gridData);


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

#endif /* HAVE_DUNE_ALUGRID */
