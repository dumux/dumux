//$Id$
#include "config.h"

#include <dune/grid/sgrid.hh>
#include "new_brineco2problem_2pni.hh"

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
        // Set the type for scalar values (should be one of float, double
        // or long double)
        typedef double                            Scalar;
        typedef Dune::NewBrineCO2Problem2pni<Scalar>  Problem;
        typedef Problem::DomainTraits::Grid       Grid;
        typedef Dune::GridPtr<Grid>               GridPointer;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments
        if (argc < 4)
            usage(argv[0]);

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
        //        typedef Dune::UGGrid<dim>                     Grid;

        // load the grid from file
        GridPointer gridPtr =  GridPointer(dgfFileName);
        Dune::gridinfo(*gridPtr);

        // instantiate and run the concrete problem
        Problem problem(&(*gridPtr), dt, tEnd);
        if (!problem.simulate())
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
