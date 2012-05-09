// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Alexandru Tatomir                            *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
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
/*!
 * \file
 *
 * \brief test for the two-phase box model
 */
#include "config.h"
#include "2pDFMMtestProblem.hh"
#include <dumux/common/start.hh>

#include <iostream>

/*!
 * \brief Provides an interface for customizing error messages associated with
 *        reading in parameters.
 *
 * \param progName  The name of the program, that was tried to be started.
 * \param errorMsg  The error message that was issued by the start function.
 *                  Comprises the thing that went wrong and a general help message.
 */
void usage(const char *progName, const std::string &errorMsg)
{
    if (errorMsg.size() > 0) {
        std::string errorMessageOut = "\nUsage: ";
                    errorMessageOut += progName;
                    errorMessageOut += " [options]\n";
                    errorMessageOut += errorMsg;
                    errorMessageOut += "\n\nThe List of Mandatory arguments for this program is:\n"
                                        "\t-tEnd                               The end of the simulation [s] \n"
                                        "\t-dtInitial                          The initial timestep size [s] \n"
                                        "\t-gridFile                           The file name of the file containing the grid \n"
                                        "\t                                        definition in DGF format\n"
                                        "\t-SpatialParameters.lensLowerLeftX   Dimension of the lens [m] \n"
                                        "\t-SpatialParameters.lensLowerLeftY   Dimension of the lens [m] \n"
                                        "\t-SpatialParameters.lensUpperRightX  Dimension of the lens [m] \n"
                                        "\t-SpatialParameters.lensUpperRighty  Dimension of the lens [m] \n"
                                        "\n";

        std::cout << errorMessageOut
                  << "\n";
    }
}

void usage(const char *progName)
{
	std::string errorMessageOut = "\nUsage: ";
				errorMessageOut += progName;
				errorMessageOut += " [--restart restartTime] grid tEnd dt\n";
    std::cout << errorMessageOut <<"\n";
    exit(1);
}

////////////////////////
// the main function
////////////////////////
int main(int argc, char** argv)
{
//    typedef TTAG(TwoPDFMMProblem) TypeTag;
//    return Dumux::start<TypeTag>(argc, argv, usage);

    try {
        FILE *file;
        typedef TTAG(TwoPDFMMTestProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
        typedef GET_PROP_TYPE(TypeTag, Grid) Grid;
        typedef GET_PROP_TYPE(TypeTag, Problem) Problem;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;

        static const int dim = Grid::dimension;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc < 4)
            usage(argv[0]);

        // deal with the restart stuff
        int argPos = 1;
        bool restart = false;
        double restartTime = 0;
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;

            std::istringstream(argv[argPos++]) >> restartTime;
        }

        if (argc - argPos != 3) {
            usage(argv[0]);
        }

        // read the initial time step and the end time
        double tEnd, dt;
        const char *artFileName = argv[argPos++];
        std::istringstream(argv[argPos++]) >> tEnd;
        std::istringstream(argv[argPos++]) >> dt;
        file = fopen(artFileName,"r");
        if (file == NULL) {
            std::cout<<"Could NOT open file"<<std::endl;
            exit(1);
        }

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////

//        Dumux::OverlapARTGeometry pT; // overlaps a grid over the fracture art file
//        pT.overlap_art_geometry();

//        Dumux::ArtReader<Grid> myARTgeometry;
//        Grid *gridT;
//        myARTgeometry.art_readARTfile(file);
//        gridT = myARTgeometry.createGrid();
//        Dumux::FractureMapper<Grid> ModelARTReader(*gridT);
//        ModelARTReader.fractureMapper(*gridT);



        ///***********
        Dumux::ArtReader<Grid> myARTgeometry;
        Grid *gridT;
        myARTgeometry.read_art_file(artFileName);
        gridT = myARTgeometry.createGrid();
        gridT->loadBalance();
//        Dumux::FractureMapper<Grid> ModelARTReader(*gridT);
        Dumux::FractureMapper<Grid> ModelARTReader(*gridT, myARTgeometry);
        ModelARTReader.fractureMapper(*gridT);
//        myARTgeometry.new_fractureMapper(*gridT);
//
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
//        Problem problem(timeManager, gridT->leafView(), myARTgeometry.isDuneFractureVertex_,myARTgeometry.isDuneFractureEdge_);
        timeManager.init(problem, 0, dt, tEnd, restart);
        if (restart)
            problem.restart(restartTime);
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

    return 3;

}
