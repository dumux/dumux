//$Id$
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Melanie Darcis                               *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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

#include "injectionproblem2pni.hh"

#include <dune/grid/common/gridinfo.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

void usage(const char *progname)
{
    std::cout << boost::format("usage: %s [--restart restartTime] tEnd dt\n")%progname;
    exit(1);
};

int main(int argc, char** argv)
{
    try {
        typedef TTAG(InjectionProblem2PNI) TypeTag;
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
        if (argc < 3)
            usage(argv[0]);

        // deal with the restart arguments
        int argPos = 1;
        bool restart = false;
        double restartTime = 0;
        // check if simulation should be initialised by a restart file
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;
            // read time for which restart should be performed
            std::istringstream(argv[argPos++]) >> restartTime;
        }

        if (argc - argPos != 2) {
            usage(argv[0]);
        }

        double tEnd, dt;
        std::istringstream(argv[argPos++]) >> tEnd;
        std::istringstream(argv[argPos++]) >> dt;

        // create the grid
        GlobalPosition lowerLeft(0.0);
        GlobalPosition upperRight;
        Dune::FieldVector<int,dim> res; // cell resolution
        upperRight[0] = 100.0;
        upperRight[1] = 50.0;
        res[0]        = 25;
        res[1]        = 20;

#if USE_UG
        ////////////////////////////////////////////////////////////
        // Make a uniform grid
        ////////////////////////////////////////////////////////////
        Grid grid;

        Dune::GridFactory<Grid> factory(&grid);
        for (int i=0; i<=res[0]; i++) {
            for (int j=0; j<=res[1]; j++) {
                Dune::FieldVector<double,2> pos;
                pos[0] = upperRight[0]*double(i)/res[0];
                pos[1] = upperRight[1]*double(j)/res[1];
                factory.insertVertex(pos);
            }
        }

        for (int i=0; i<res[0]; i++) {
            for (int j=0; j<res[1]; j++) {
                std::vector<unsigned int> v(4);
                v[0] = i*(res[1]+1) + j;
                v[1] = i*(res[1]+1) + j+1;
                v[2] = (i+1)*(res[1]+1) + j;
                v[3] = (i+1)*(res[1]+1) + j+1;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), v);
            }
        }

        factory.createGrid();
#else
        Grid grid(
#ifdef HAVE_MPI
                  Dune::MPIHelper::getCommunicator(),
#endif
                  upperRight, // upper right
                  res, // number of cells
                  Dune::FieldVector<bool,dim>(false), // periodic
                  2); // overlap
#endif

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        // specify dimensions of the low-permeable lens
        GlobalPosition lowerLeftLens, upperRightLens;
        lowerLeftLens[0] = 0.0;
        lowerLeftLens[1] = 30.0;
        upperRightLens[0] = 100.0;
        upperRightLens[1] = 30.0;
        Problem problem(grid.leafView(),
                        lowerLeftLens, upperRightLens);
        
        // load restart file if necessary
        if (restart)
            problem.deserialize(restartTime);
        // run the simulation
        if (!problem.simulate(dt, tEnd))
            return 2;

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
