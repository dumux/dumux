#include "config.h"

#include "new_benchmark3problem.hh"
#include "new_easyproblem.hh"
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser.hh>
#include <dune/grid/uggrid.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

void usage(const char *progname)
{
    std::cout << boost::format("usage: %s [--restart restartTime] basefilename tEnd dt\n")%progname;
    exit(1);
};

int main(int argc, char** argv)
{
    try {
        // Set the type for scalar values (should be one of float, double
        // or long double)
        typedef double                            Scalar;
        typedef Dune::NewBenchmark3Problem<Scalar>        Problem;
        typedef Problem::DomainTraits::Grid       Grid;
        typedef Dune::GridPtr<Grid>               GridPointer;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
        if (argc < 4)
            usage(argv[0]);

        double tEnd, dt;
        bool restart = false;
        int argPos = 1;
        double restartTime = 0;
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;

            std::istringstream(argv[argPos++]) >> restartTime;
        }

        if (argc - argPos != 3) {
            usage(argv[0]);
        }

        const char *dgfFileName = argv[argPos++];
        std::istringstream(argv[argPos++]) >> tEnd;
        std::istringstream(argv[argPos++]) >> dt;

        // load the grid from file
        GridPointer gridPtr =  GridPointer(dgfFileName);
        Dune::gridinfo(*gridPtr);

        // instantiate and run the concrete problem
        Problem problem(&(*gridPtr), dt, tEnd);

        // load restart file if necessarry
        if (restart)
            problem.deserialize(restartTime);

        if (!problem.simulate())
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
