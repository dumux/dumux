// $Id: test_1p2c.cc 3779 2010-06-24 07:07:56Z bernd $
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
#include "config.h"

#include "tissue_tumor_problem.hh"

#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridinfo.hh>

#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

void usage(const char *progname)
{
    std::cout << boost::format("usage: %s [--restart restartTime] gridFile.dgf tEnd dt\n")%progname;
    exit(1);
};

int main(int argc, char** argv)
{
    try {
        typedef TTAG(TissueTumorProblem)              TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar))  Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid))    Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;
        typedef Dune::GridPtr<Grid>                             GridPointer;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
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

        double tEnd, dt;
        const char *dgfFileName = argv[argPos++];
        std::istringstream(argv[argPos++]) >> tEnd;
        std::istringstream(argv[argPos++]) >> dt;

        // create grid

        // -> load the grid from file
        GridPointer gridPtr =  GridPointer(dgfFileName);
        (*gridPtr).loadBalance();
        Dune::gridinfo(*gridPtr);

        // instantiate and run the concrete problem
        Problem problem(gridPtr->leafView());
        if (!problem.simulate(dt, tEnd))
            return 2;

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
