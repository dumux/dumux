// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
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
 * \brief Provides a few default main functions for convenience.
 */
#ifndef DUMUX_START_HH
#define DUMUX_START_HH

#include <dumux/common/propertysystem.hh>

#include <dune/grid/io/file/dgfparser.hh>
#include <dune/common/mpihelper.hh>
#include <iostream>

namespace Dumux
{
// forward declaration of property tags
namespace Properties
{
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(TimeManager);
}

/*!
 * \brief Print a usage string for simulations using
 *        Dumux::startFromDGF() as their main() function.
 *
 * \param progname The name of the executable
 */
void printUsageDGF(const char *progname)
{
    std::cout << "usage: " << progname << " [--restart restartTime] gridFile.dgf tEnd dt\n";
    exit(1);
};

/*!
 * \brief Print a usage string for simulations using
 *        Dumux::startWithGrid() as their main() function.
 *
 * \param progname The name of the executable
 */
void printUsageGrid(const char *progname)
{
    std::cout << "usage: " << progname << " [--restart restartTime] tEnd dt\n";
    exit(1);
};

/*!
 * \brief Provides a default main function for simulations requiring
 *        only a single DGF file as their grid specification.
 *
 * \tparam ProblemTypeTag  The type tag of the problem which needs to be solved
 *
 * \param argc  The 'argc' argument of the main function
 * \param argv  The 'argv' argument of the main function
 */
template <class ProblemTypeTag>
int startFromDGF(int argc, char **argv)
{
#ifdef NDEBUG
    try {
#endif

        typedef typename GET_PROP_TYPE(ProblemTypeTag, PTAG(Grid)) Grid;
        typedef typename GET_PROP_TYPE(ProblemTypeTag, PTAG(Problem)) Problem;
        typedef typename GET_PROP_TYPE(ProblemTypeTag, PTAG(TimeManager)) TimeManager;
        typedef Dune::GridPtr<Grid> GridPointer;

        // initialize MPI, finalize is done automatically on exit
        static Dune::MPIHelper& mpiHelper = Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
        if (argc < 4)
            printUsageDGF(argv[0]);

        // deal with the restart stuff
        int argIdx = 1;
        bool restart = false;
        double restartTime = 0;
        if (std::string("--restart") == argv[argIdx]) {
            restart = true;
            ++argIdx;

            std::istringstream(argv[argIdx++]) >> restartTime;
        }

        if (argc - argIdx != 3) {
            printUsageDGF(argv[0]);
        }

        double tEnd, dt;
        const char *dgfFileName = argv[argIdx++];
        std::istringstream(argv[argIdx++]) >> tEnd;
        std::istringstream(argv[argIdx++]) >> dt;

        // create grid
        // -> load the grid from file
        GridPointer gridPtr(dgfFileName);
        if (mpiHelper.size() > 1) {
            if (!Dune::Capabilities::isParallel<Grid>::v) {
                std::cerr << "WARNING: THE PROGRAM IS STARTED USING MPI, BUT THE GRID IMPLEMENTATION\n"
                          << "         YOU HAVE CHOSEN IS NOT PARALLEL!\n";
            }
            gridPtr.loadBalance();
        }

        // instantiate and run the concrete problem
        TimeManager timeManager;
        Problem problem(timeManager, gridPtr->leafView());
        timeManager.init(problem, 0, dt, tEnd, !restart);
        if (restart)
            problem.restart(restartTime);
        timeManager.run();
        return 0;

#ifdef NDEBUG
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
    }
#endif

    return 3;
};


/*!
 * \brief Provides a default main function for simulations which
 *        create the grid themself but do not require any other
 *        parameters.
 *
 * \tparam ProblemTypeTag  The type tag of the problem which needs to be solved
 *
 * \param grid  The grid used by the simulation
 * \param argc  The 'argc' argument of the main function
 * \param argv  The 'argv' argument of the main function
 */
template <class TypeTag>
int startWithGrid(const typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) &grid,
                  int argc,
                  char **argv)
{
#ifdef NDEBUG
    try {
#endif

        typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef typename GET_PROP_TYPE(TypeTag, PTAG(TimeManager)) TimeManager;

        // parse the command line arguments for the program
        if (argc < 3)
            printUsageGrid(argv[0]);

        // deal with the restart stuff
        int argIdx = 1;
        bool restart = false;
        double restartTime = 0;
        if (std::string("--restart") == argv[argIdx]) {
            restart = true;
            ++argIdx;

            std::istringstream(argv[argIdx++]) >> restartTime;
        }

        if (argc - argIdx != 2) {
            printUsageGrid(argv[0]);
        }

        double tEnd, dt;
        std::istringstream(argv[argIdx++]) >> tEnd;
        std::istringstream(argv[argIdx++]) >> dt;

        // instantiate and run the concrete problem
        TimeManager timeManager;
        Problem problem(timeManager, grid.leafView());
        timeManager.init(problem, 0, dt, tEnd, !restart);
        if (restart)
            problem.restart(restartTime);
        timeManager.run();
        return 0;

#ifdef NDEBUG
    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
    }
#endif

    return 3;
};

}

#endif
