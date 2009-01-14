#include "config.h"
#include "new_waterairproblem.hh"

#include <dune/common/exceptions.hh>

#include <iostream>

int main(int argc, char** argv)
{
    // Set the type for scalar values (should be one of float, double
    // or long double)
    typedef double                            Scalar;
    typedef Dune::NewWaterAirProblem<Scalar>  Problem;
    typedef Problem::DomainTraits::Grid       Grid;
    typedef Dune::GridPtr<Grid>               GridPointer;

    try {
        // parse the command line arguments for the program
        if (argc != 4) {
            std::cout << boost::format("usage: %s gridFile.dgf tEnd dt\n")%argv[0];
            return 1;
        }
        double tEnd, dt;
        const char *dgfFileName = argv[1];
        std::istringstream(argv[2]) >> tEnd;
        std::istringstream(argv[3]) >> dt;

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
    }
    return 3;
}
