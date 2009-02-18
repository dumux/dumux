#include "config.h"

#include "new_1pproblem.hh"

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfalberta.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

int main(int argc, char** argv)
{
    try {
        // Set the type for scalar values (should be one of float, double
        // or long double)
        const int dim = 3;
        typedef double                                   Scalar;
        typedef Dune::ALUCubeGrid<dim, dim>              Grid;
//        typedef Dune::YaspGrid<dim>                      Grid;
//        typedef Dune::UGGrid<dim>                        Grid;
        typedef Dune::NewOnePProblem<Grid, Scalar>       Problem;
        typedef Problem::DomainTraits::GlobalPosition    GlobalPosition;
        typedef Dune::GridPtr<Grid>                      GridPointer;
        
        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
        if (argc != 4) {
            std::cout << boost::format("usage: %s dgfFileName tEnd dt\n")%argv[0];
            return 1;
        }
        double tEnd, dt;
        const char *dgfFileName = argv[1];
        std::istringstream(argv[2]) >> tEnd;
        std::istringstream(argv[3]) >> dt;

        // create the grid from a DGF file
        GridPointer gridPtr = GridPointer(dgfFileName,
                                          Dune::MPIHelper::getCommunicator());
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
