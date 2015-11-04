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
#include "config.h"
#include <iostream>
#include <dune/common/parametertreeparser.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dumux/io/cakegridcreator.hh>
#include <dumux/common/basicproperties.hh>

namespace Dumux
{

namespace Properties
{
NEW_TYPE_TAG(GridCreatorCakeTest, INHERITS_FROM(NumericModel));
//     Set the grid type
#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(GridCreatorCakeTest, Grid, Dune::ALUGrid<3, 3, Dune::cube, Dune::conforming>);
#elif HAVE_UG
SET_TYPE_PROP(GridCreatorCakeTest, Grid, Dune::UGGrid<3>);
#else
SET_TYPE_PROP(GridCreatorCakeTest, Grid, Dune::YaspGrid<3>);
#endif
}
}

int main(int argc, char** argv)
{
#if HAVE_DUNE_ALUGRID
    try {
        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // Some typedefs
        typedef typename TTAG(GridCreatorCakeTest) TypeTag;
        typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef typename Dumux::LinearSpacer<TypeTag> LinearSpacer;
        typedef typename Dumux::CakeGridCreator<TypeTag> GridCreatorLinear;
        typedef typename Dumux::PowerLawSpacer<TypeTag> PowerLawSpacer;
        typedef typename Dumux::CakeGridCreator<TypeTag, PowerLawSpacer> GridCreatorPowerLaw;
        //    typedef typename Dumux::GridCreator<TypeTag> GridCreator;

        // Read the parameters from the input file
        typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
        Dune::ParameterTreeParser::readINITree("test_gridcreator_cake.input", ParameterTree::tree());

//      Make the grid
        GridCreatorLinear::makeGrid();
        GridCreatorPowerLaw::makeGrid();

        // construct a vtk output writer and attach the boundaryMakers
        Dune::VTKSequenceWriter<Grid::LeafGridView> vtkWriterLinear(GridCreatorLinear::grid().leafGridView(), "pieceofcake-linspace", ".", "");
        //    vtkWriter.addVertexData(boundaryMarker, "boundaryMarker");
        vtkWriterLinear.write(0);

        Dune::VTKSequenceWriter<Grid::LeafGridView> vtkWriterPowerLaw(GridCreatorPowerLaw::grid().leafGridView(), "cake-powerspace", ".", "");
        //    vtkWriter.addVertexData(boundaryMarker, "boundaryMarker");
        vtkWriterPowerLaw.write(0);

        return 0;
    }
    catch (Dumux::ParameterException &e) {
        typedef typename TTAG(GridCreatorCakeTest) TypeTag;
        Dumux::Parameters::print<TypeTag>();
        std::cerr << e << ". Abort!\n";
        return 1;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
        return 3;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        return 4;
    }
#else
#warning "You need to have ALUGrid or UGGrid installed to run this test."
    std::cerr << "You need to have ALUGrid or UGGrid installed to run this test\n";
    return 77;
#endif
} //closing main
//} //closing namespace dumux

