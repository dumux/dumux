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
 * \brief test for the sequential 2p2c model
 */
#include "config.h"

#include "test_dec2p2cproblem.hh"

#include <dune/grid/common/gridinfo.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <iostream>

////////////////////////
// the main function
////////////////////////
void usage(const char *progname)
{
    std::cout << "usage: " << progname << " [--restart restartTime] tEnd firstDt\n";
    exit(1);
}

int main(int argc, char** argv)
{
    try {
        typedef TTAG(TestDecTwoPTwoCProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar))  Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid))    Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;

        static const int dim = Grid::dimension;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        // deal with the restart stuff
        int argPos = 1;
        bool restart = false;
        double startTime = 0;
        // deal with start parameters
        double tEnd= 3e3;
        double firstDt = 200;
        if (argc != 1)
        {
            // deal with the restart stuff
            if (std::string("--restart") == argv[argPos]) {
                restart = true;
                ++argPos;

                std::istringstream(argv[argPos++]) >> startTime;
            }
            if (argc - argPos == 2)
            {
                // read the initial time step and the end time
                std::istringstream(argv[argPos++]) >> tEnd;
                std::istringstream(argv[argPos++]) >> firstDt;
            }
            else
                usage(argv[0]);
        }
        else
        {
            Dune::dwarn << "simulation started with predefs" << std::endl;
        }

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////
        Dune::FieldVector<int,dim> N(10);
        Dune::FieldVector<double ,dim> L(0);
        Dune::FieldVector<double,dim> H(10);
        Grid grid(
#ifdef HAVE_MPI
            Dune::MPIHelper::getCommunicator(),
#endif
            H, // upper right
            N, // number of cells
            Dune::FieldVector<bool,dim>(false), // periodic
            1); // overlap

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        TimeManager timeManager;
        Problem problem(timeManager, grid.leafView(), L, H);

        // initialize the simulation
        timeManager.init(problem, startTime, firstDt, tEnd, restart);
        // run the simulation
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
