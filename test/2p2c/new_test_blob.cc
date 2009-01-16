#define ENABLE_MPI 1

#include "config.h"
#include "new_blobproblem.hh"

#include <dune/common/exceptions.hh>

#include <iostream>

int main(int argc, char** argv)
{
    // Set the type for scalar values (should be one of float, double
    // or long double)
    typedef double                            Scalar;
    typedef Dune::NewBlobProblem<Scalar>      Problem;
    typedef Problem::DomainTraits::Grid       Grid;
    typedef Problem::DomainTraits::GlobalPosition GlobalPosition;
    typedef Dune::GridPtr<Grid>               GridPointer;

    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    for (int i = 0; i < argc; ++i)
        std::cout << argv[i] << "\n";
    
    try {
        // parse the command line arguments for the program
        if (argc != 3) {
            std::cout << boost::format("usage: %s tEnd dt\n")%argv[0];
            return 1;
        }

        double tEnd, dt;
        std::istringstream(argv[1]) >> tEnd;
        std::istringstream(argv[2]) >> dt;

        // instantiate and run the concrete problem
        Dune::NewBlobProblem<Scalar> problem(dt, tEnd);
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
