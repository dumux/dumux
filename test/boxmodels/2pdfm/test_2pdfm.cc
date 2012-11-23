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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Test for the two-phase flow with discrete fracture-matrix (discrete
 *        fracture model) and box model scheme.
 */
#include "config.h"

#include "2pdfmtestproblem.hh"
#include <dumux/common/start.hh>

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 */
void usage(const char *progName)
{
    std::string errorMessageOut = "\nUsage: ";
                errorMessageOut += progName;
                errorMessageOut += " [--restart restartTime] grid tEnd dt\n";
    std::cout << errorMessageOut << std::endl;
    exit(1);
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
    try {
        FILE *file;
        typedef TTAG(TwoPDFMTestProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
        typedef GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc < 4)
        {
            usage(argv[0]);
        }

        // deal with the restart stuff
        int argPos = 1;
        bool restart = false;
        double restartTime = 0;
        if (std::string("--restart") == argv[argPos])
        {
            restart = true;
            ++argPos;

            std::istringstream(argv[argPos++]) >> restartTime;
        }

        if (argc - argPos != 3)
        {
            usage(argv[0]);
        }

        // read the initial time step and the end time
        double tEnd, dt;
        const char *artFileName = argv[argPos++];
        std::istringstream(argv[argPos++]) >> tEnd;
        std::istringstream(argv[argPos++]) >> dt;
        file = fopen(artFileName,"r");
        if (file == NULL)
        {
            std::cout << "Could not open artmesh file '" 
                      << artFileName << "'" << std::endl;
            exit(1);
        }

        ////////////////////////////////////////////////////////////
        // create the grid from ART Reader
        ////////////////////////////////////////////////////////////

        Dumux::ArtReader<Grid> dfmArtmeshGeometry; // instantiate ArtReader
        Grid *gridT;
        dfmArtmeshGeometry.read_art_file(artFileName);
        gridT = dfmArtmeshGeometry.createGrid();
        gridT->loadBalance();
        Dumux::FractureMapper<Grid> ModelARTReader(*gridT, dfmArtmeshGeometry); // map fractures to grid
        ModelARTReader.fractureMapper();

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        TimeManager timeManager;
        Problem problem(timeManager,
                            gridT->leafView(),
                            ModelARTReader.isDuneFractureVertex_,
                            ModelARTReader.isDuneFractureEdge_,
                            ModelARTReader.fractureEdgesIdx_);
        timeManager.init(problem, 0, dt, tEnd, restart);
        if (restart)
        {
            problem.restart(restartTime);
        }
        timeManager.run();
        return 0;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }
    return EXIT_SUCCESS;
}
