#include "config.h"

#include <dune/grid/sgrid.hh>
#include "seqcoup_base_problem.hh"
#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

void usage(const char *progname)
{
    std::cout << boost::format("usage: %s [--restart restartTime] grid tEnd_FirstModel dt_FirstModel tEnd_SecondModel dt_SecondModel\n")%progname;
    exit(1);
};

int main(int argc, char** argv)
{
    try {
        // Set the type for scalar values (should be one of float, double
        // or long double)
        typedef double                                  Scalar;
//        typedef Dune::SeqCoup2P2CNIProblem<Scalar>      Problem2P2CNI;
//        typedef Dune::SeqCoup2PNIProblem<Scalar>        Problem2PNI;
        typedef Dune::SeqCoupBaseProblem<Scalar>        ProblemBase;
        typedef ProblemBase::DomainTraits::Grid         Grid;
        typedef Dune::GridPtr<Grid>                     GridPointer;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments
        if (argc < 6)
            usage(argv[0]);

        int argPos = 1;
        bool restart = false;
        double restartTime = 0;
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;

            std::istringstream(argv[argPos++]) >> restartTime;
        }

        if (argc - argPos != 5) {
            usage(argv[0]);
        }
        double tEnd2pni,tEnd2p2cni, dt2pni, dt2p2cni;
        const char *dgfFileName = argv[argPos++];
        std::istringstream(argv[argPos++]) >> tEnd2pni;
        std::istringstream(argv[argPos++]) >> dt2pni;
        std::istringstream(argv[argPos++]) >> tEnd2p2cni;
        std::istringstream(argv[argPos++]) >> dt2p2cni;

        // load the grid from file
        GridPointer gridPtr =  GridPointer(dgfFileName);
        Dune::gridinfo(*gridPtr);

        // instantiate and run the concrete problem
        ProblemBase problembase(&(*gridPtr), dt2pni, tEnd2pni,dt2p2cni, tEnd2p2cni);

        // load restart file if necessarry
        if (restart)
             problembase.deserialize(restartTime);
        if (!problembase.simulate())
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
