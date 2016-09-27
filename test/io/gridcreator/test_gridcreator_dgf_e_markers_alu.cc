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
 * \brief Test for gmsh interface of the grid creator
 */
#include <config.h>
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dumux/common/basicproperties.hh>

#include "gridcreatortests.hh"

namespace Dumux {
namespace Properties {

    NEW_TYPE_TAG(GridCreatorGmshTest, INHERITS_FROM(NumericModel));
#if HAVE_DUNE_ALUGRID
    SET_TYPE_PROP(GridCreatorGmshTest, Grid, Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>);
    SET_BOOL_PROP(GridCreatorGmshTest, AdaptiveGrid, false);
#else
    SET_TYPE_PROP(GridCreatorGmshTest, Grid, Dune::YaspGrid<2>);
#endif
    // Change the default "Grid" to customized "BifurcationGrid", merely for demonstration purposes.
    SET_STRING_PROP(GridCreatorGmshTest, GridParameterGroup, "BifurcationGrid");
}
}

int main(int argc, char** argv)
{
#if HAVE_DUNE_ALUGRID
try {
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);
    Dumux::GridCreatorTests<TTAG(GridCreatorGmshTest)>::testElementDomainMarkers("test_gridcreator_dgf_e_markers.input", "dgf", "co2_alu");
    return 0;
}
catch (Dumux::ParameterException &e) {
    Dumux::Parameters::print<TTAG(GridCreatorGmshTest)>();
    std::cerr << e << ". Abort!\n";
    return 1;
}
catch (Dune::Exception &e) {
    std::cerr << "Dune reported error: " << e << std::endl;
    return 3;
}
#else
#warning "You need to have dune-alugrid installed to run this test."
    std::cerr << "You need to have dune-alugrid installed to run this test\n";
    return 77;
#endif
}
