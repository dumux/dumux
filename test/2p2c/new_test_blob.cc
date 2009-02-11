#include "config.h"

#include "new_blobproblem.hh"

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser.hh>
#include <dune/grid/uggrid.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

int main(int argc, char** argv)
{
    try {

        // Set the type for scalar values (should be one of float, double
        // or long double)
        const int dim = 2;
        typedef double                                Scalar;
//        typedef Dune::YaspGrid<dim>                   Grid;
//        typedef Dune::ALUSimplexGrid<dim,dim>         Grid;
        typedef Dune::UGGrid<dim>                     Grid;
        typedef Dune::GridPtr<Grid>                   GridPointer;
        typedef Dune::NewBlobProblem<Grid, Scalar>    Problem;
        typedef Problem::DomainTraits::GlobalPosition GlobalPosition;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
        if (argc != 4) {
            std::cout << boost::format("usage: %s tEnd dt dgfFile\n")%argv[0];
            return 1;
        }

        double tEnd, dt;
        std::istringstream(argv[1]) >> tEnd;
        std::istringstream(argv[2]) >> dt;

        const char *dgfFileName = argv[3];
        
        // load the grid from file
        GridPointer gridPtr = GridPointer(dgfFileName,
                                          Dune::MPIHelper::getCommunicator());
        Dune::gridinfo(*gridPtr);

/*
        Grid grid(
#ifdef HAVE_MPI
                  Dune::MPIHelper::getCommunicator(),
#endif
                  GlobalPosition(300.0), // upper right
                  Dune::FieldVector<int,dim>(40), // number of cells
                  Dune::FieldVector<bool,dim>(false), // periodic
                  2); // overlap
*/

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
