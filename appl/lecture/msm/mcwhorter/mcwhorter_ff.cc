/*****************************************************************************
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#include "config.h"
#include <iostream>
#include <boost/format.hpp>
#include <iomanip>

#include <dune/grid/common/gridinfo.hh>
#include <dune/common/parametertreeparser.hh>
#include "mcwhorterproblem.hh"

void usage(const char *progname)
{
    std::cout << boost::format("usage: %s InputFileName\n")%progname;
    exit(1);
}

int main(int argc, char** argv)
{
    try
    {
        typedef TTAG(McWhorterProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;
        typedef GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

        static const int dim = Grid::dimension;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc != 2)
            usage(argv[0]);

        std::string inputFileName;
        inputFileName = argv[1];

        ////////////////////////////////////////////////////////////
        // Read Input file
        ////////////////////////////////////////////////////////////

        Dune::ParameterTree inputParameters;
        Dune::ParameterTreeParser::readINITree(inputFileName, Params::tree());

        double discretizationLength = Params::tree().get<double>("problem.DiscretizationLength");

         int cellNumber = static_cast<int>(2/discretizationLength);

        // define the problem dimensions
        Dune::FieldVector<Scalar, dim> lowerLeft(0);
        Dune::FieldVector<Scalar, dim> upperRight(2.0);
//        UpperRight[1] = 1;
        Dune::FieldVector<int, dim> cellNumbers(cellNumber);
//        cellNumbers[0] = 26;

        double tEnd = Params::tree().get<double>("problem.tEnd");
        double dt = tEnd;

        // create a grid object
        typedef Dune::SGrid<dim, dim> GridType;
        typedef GridType::LeafGridView GridView;

        // grid reference
        GridType grid(cellNumbers,lowerLeft,upperRight);
//        Dune::gridinfo(grid);


        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////

        typedef GET_PROP_TYPE(TypeTag, PTAG(WettingPhase)) WettingPhase;
        typedef GET_PROP_TYPE(TypeTag, PTAG(NonWettingPhase)) NonWettingPhase;

        WettingPhase::Component::setViscosity(Params::tree().get<double>("fluid.viscosityW"));
        NonWettingPhase::Component::setViscosity(Params::tree().get<double>("fluid.viscosityNW"));

        TimeManager timeManager;
        Problem problem(timeManager, grid.leafView(), lowerLeft, upperRight);

        timeManager.init(problem, 0, dt, tEnd);
        timeManager.run();

        return 0;

    } catch (Dune::Exception &e)
    {
        std::cerr << "Dune reported error: " << e << std::endl;
    } catch (...)
    {
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
