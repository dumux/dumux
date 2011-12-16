// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 20010 by Markus Wolff                                     *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
 * \ingroup IMPETtests
 * \brief test for the sequential 2p model
 */
#include "config.h"

#include "test_impes_problem.hh"

#include <dune/grid/common/gridinfo.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>


////////////////////////
// the main function
////////////////////////
void usage(const char *progname)
{
    std::cout << boost::format("usage: %s [--restart restartTime] tEnd\n")%progname;
    exit(1);
}

int main(int argc, char** argv)
{
    try {
        typedef TTAG(IMPESTestProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;

        static const int dim = Grid::dimension;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc < 2)
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

        if (argc - argPos != 1) {
            usage(argv[0]);
        }

        // read the initial time step and the end time
        double tEnd, dt;
        std::istringstream(argv[argPos++]) >> tEnd;
        dt = tEnd;

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////
        Dune::FieldVector<int,dim> numCells;
        numCells[0] = 30;
        numCells[1] = 6;
        Dune::FieldVector<double,dim> lowerLeft;
        lowerLeft[0] = 0;
        lowerLeft[1] = 0;
        Dune::FieldVector<double,dim> upperRight;
        upperRight[0] = 300;
        upperRight[1] = 60;
        Dune::FieldVector<bool,dim> periodic(false);

        Grid grid(Dune::MPIHelper::getCommunicator(), upperRight, numCells, periodic, /*overlap=*/1);
        grid.loadBalance();
        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        TimeManager timeManager;
        Problem problem(timeManager, grid.leafView());

        // use restart file if necessarry
        timeManager.init(problem, restartTime, dt, tEnd, restart);
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
