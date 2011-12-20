// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
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
 * \brief test for the Richards box model
 */
#include "config.h"

#include "richardslensproblem.hh"
#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridinfo.hh>

#include <dune/common/mpihelper.hh>

#include <iostream>

void usage(const char *progname)
{
    std::cout << "usage: " << progname << " [--restart restartTime] gridFile.dgf tEnd dt\n";
    exit(1);
};

int main(int argc, char** argv)
{
    try {
        typedef TTAG(RichardsLensProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;
        typedef Dune::GridPtr<Grid> GridPointer;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
        if (argc < 4)
            usage(argv[0]);

        // deal with the restart stuff
        int argPos = 1;
        bool restart = false;
        double startTime = 0;
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;

            std::istringstream(argv[argPos++]) >> startTime;
        }

        if (argc - argPos != 3) {
            usage(argv[0]);
        }

        double tEnd, dt;
        const char *dgfFileName = argv[argPos++];
        std::istringstream(argv[argPos++]) >> tEnd;
        std::istringstream(argv[argPos++]) >> dt;

        // create grid

        // -> load the grid from file
        GridPointer gridPtr = GridPointer(dgfFileName);
        (*gridPtr).loadBalance();
       Dune::gridinfo(*gridPtr);


        // specify dimensions of the low-permeable lens
        GlobalPosition lowerLeftLens, upperRightLens;
        lowerLeftLens[0] = 1.0;
        lowerLeftLens[1] = 2.0;
        upperRightLens[0] = 4.0;
        upperRightLens[1] = 3.0;

        // instantiate and run the concrete problem
        TimeManager timeManager;
        Problem problem(timeManager, gridPtr->leafView(), lowerLeftLens, upperRightLens);
        timeManager.init(problem, startTime, dt, tEnd, restart);
        timeManager.run();
        return 0;
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
    }

    return 3;
}
