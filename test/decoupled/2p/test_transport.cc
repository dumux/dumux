// $Id: test_2p.cc 3732 2010-06-11 13:27:20Z bernd $
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
 * \brief test for the explicit transport model
 */
#include "config.h"

#include "test_transport_problem.hh"

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
    std::cout << boost::format("usage: %s [--restart restartTime] gridFile.dgf tEnd\n")%progname;
    exit(1);
}

int main(int argc, char** argv)
{
    try {
        typedef TTAG(TransportTestProblem) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;

        static const int dim = Grid::dimension;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        ////////////////////////////////////////////////////////////
        // parse the command line arguments
        ////////////////////////////////////////////////////////////
        if (argc < 3)
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

        if (argc - argPos != 2) {
            usage(argv[0]);
        }

        ////////////////////////////////////////////////////////////
        // create the grid
        ////////////////////////////////////////////////////////////
        const char *dgfFileName = argv[argPos++];
        Dune::GridPtr<Grid> gridPtr(dgfFileName);

        // read the initial time step and the end time
        double tEnd, dt;
        std::istringstream(argv[argPos++]) >> tEnd;
        dt = tEnd;

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        Dune::FieldVector<double ,dim> L(0);
        Dune::FieldVector<double,dim> H(1);
        Problem problem(gridPtr->leafView(), L, H);

        // load restart file if necessarry
        if (restart)
            problem.deserialize(restartTime);

        problem.timeManager().init(problem, 0, dt, tEnd, !restart);
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
