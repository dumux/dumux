// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup IMPETtests
 * \brief test for the multiphysics 2p2c model
 */
#include <config.h>

#include <array>
#include <iostream>

#include <dumux/common/properties.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/common/gridinfo.hh>

#include "test_multiphysics2p2cproblem.hh"

////////////////////////
// the main function
////////////////////////
void usage(const char *progname, const std::string &errorMsg = "")
{
    std::cout << "usage: " << progname << " [--restart restartTime] tEnd firstDt\n";
    exit(1);
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    using TypeTag = Properties::TTag::TestMultTwoPTwoC;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using TimeManager = GetPropType<TypeTag, Properties::TimeManager>;

    static const int dim = Grid::dimension;

    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    auto defaultParams = [] (Dune::ParameterTree& p) {GetProp<TypeTag, Properties::ModelDefaultParameters>::defaultParams(p);};
    Dumux::Parameters::init(argc, argv, defaultParams, usage);

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
    std::array<int,dim> N;
    std::fill(N.begin(), N.end(), 10);
    Dune::FieldVector<double,dim> H(10.0);
    Grid grid(H, N);

    ////////////////////////////////////////////////////////////
    // instantiate and run the concrete problem
    ////////////////////////////////////////////////////////////
    TimeManager timeManager;
    Problem problem(timeManager, grid, H);

    // initialize the simulation
    timeManager.init(problem, startTime, firstDt, tEnd, restart);
    // run the simulation
    timeManager.run();
    return 0;
}
