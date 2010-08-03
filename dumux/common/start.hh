// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief Provides a default main function for simulations.
 */
#ifndef DUMUX_START_HH
#define DUMUX_START_HH

#include <dumux/common/propertysystem.hh>

#include <dune/common/mpihelper.hh>

#include <dune/grid/common/gridinfo.hh>

#include <boost/format.hpp>

#include <iostream>

namespace Dumux
{
void printUsage(const char *progname)
{
    std::cout << boost::format("usage: %s [--restart restartTime] gridFile.dgf tEnd dt\n")%progname;
    exit(1);
};

/*!
 * \brief Provides a default main function for simulations.
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

        typedef typename GET_PROP_TYPE(ProblemTypeTag, PTAG(Scalar))  Scalar;
        typedef typename GET_PROP_TYPE(ProblemTypeTag, PTAG(Grid))    Grid;
        typedef typename GET_PROP_TYPE(ProblemTypeTag, PTAG(Problem)) Problem;
        typedef typename GET_PROP_TYPE(ProblemTypeTag, PTAG(TimeManager)) TimeManager;
        typedef Dune::GridPtr<Grid>                                   GridPointer;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
        if (argc < 4)
            printUsage(argv[0]);

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
            printUsage(argv[0]);
        }

        double tEnd, dt;
        const char *dgfFileName = argv[argIdx++];
        std::istringstream(argv[argIdx++]) >> tEnd;
        std::istringstream(argv[argIdx++]) >> dt;

        // create grid
        // -> load the grid from file
        GridPointer gridPtr =  GridPointer(dgfFileName);
        Dune::gridinfo(*gridPtr);

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
}

#endif
