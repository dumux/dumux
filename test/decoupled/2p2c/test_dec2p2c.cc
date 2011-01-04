// $Id: test_2pinjection.cc 3151 2010-02-15 11:41:23Z mwolff $
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
#include <boost/format.hpp>

////////////////////////
// the main function
////////////////////////
void usage(const char *progname)
{
    std::cout << boost::format("usage: %s [--restart restartTime] tEnd firstDt\n")%progname;
    exit(1);
}

int main(int argc, char** argv)
{
    try {
        typedef TTAG(TestDecTwoPTwoCProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar))  Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid))    Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
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
        double restartTime = 0;
        // deal with start parameters
        double tEnd= 3e3;
        double firstDt = 1e3;
        if (argc != 1)
        {
            // deal with the restart stuff
            if (std::string("--restart") == argv[argPos]) {
                restart = true;
                ++argPos;

                std::istringstream(argv[argPos++]) >> restartTime;
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
        Grid grid(N,L,H);

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////

        Problem problem(grid.leafView(), L, H);

        // load restart file if necessarry
        if (restart)
            problem.deserialize(restartTime);

        // initialize the simulation
        problem.timeManager().init(problem, 0, firstDt, tEnd, !restart);
        // run the simulation
        problem.timeManager().run();
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
