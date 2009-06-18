// $Id:$
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf                                     *
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

#include "injectionproblem.hh"

#include <dune/common/exceptions.hh>
#include <dune/grid/common/gridinfo.hh>

#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

int main(int argc, char** argv)
{
    try {
        typedef GET_PROP_TYPE(TTAG(InjectionProblem), PTAG(Scalar))  Scalar;
        typedef GET_PROP_TYPE(TTAG(InjectionProblem), PTAG(Grid))    Grid;
        typedef GET_PROP_TYPE(TTAG(InjectionProblem), PTAG(Problem)) Problem;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;
        typedef Dune::GridPtr<Grid>                             GridPointer;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
        if (argc != 4) {
            std::cout << boost::format("usage: %s gridFile.dgf tEnd dt\n")%argv[0];
            return 1;
        }
        double tEnd, dt;
        const char *dgfFileName = argv[1];
        std::istringstream(argv[2]) >> tEnd;
        std::istringstream(argv[3]) >> dt;

        // create grid

        // -> load the grid from file
        GridPointer gridPtr =  GridPointer(dgfFileName);
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
