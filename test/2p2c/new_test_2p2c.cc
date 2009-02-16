#include "config.h"

#include "new_injectionproblem.hh"

#include <dune/grid/common/gridinfo.hh>
//#include <dune/grid/io/file/dgfparser.hh>

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
        const int dim = 2;
        typedef double                                   Scalar;
//        typedef Dune::ALUSimplexGrid<dim, dim>           Grid;
        typedef Dune::YaspGrid<dim>                      Grid;
//        typedef Dune::UGGrid<dim>                        Grid;
        typedef Dune::NewInjectionProblem<Grid, Scalar>  Problem;
        typedef Problem::DomainTraits::GlobalPosition    GlobalPosition;
//        typedef Dune::GridPtr<Grid>                      GridPointer;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
        if (argc != 3) {
            std::cout << boost::format("usage: %s tEnd dt\n")%argv[0];
            return 1;
        }
        double tEnd, dt;
        std::istringstream(argv[1]) >> tEnd;
        std::istringstream(argv[2]) >> dt;

        // create grid
        GlobalPosition upperRight;
        Dune::FieldVector<int,dim> res; // cell resolution
        upperRight[0] = 60.0;
        res[0]        = 24;

        upperRight[1] = 40.0;
        res[1]        = 16;

        Grid grid(
//#ifdef HAVE_MPI
//                  Dune::MPIHelper::getCommunicator(),
//#endif
                  upperRight, // upper right
                  res, // number of cells
                  Dune::FieldVector<bool,dim>(false), // periodic
                  2); // overlap
/*

        // load the grid from file
        GridPointer gridPtr =  GridPointer(dgfFileName,
                                           Dune::MPIHelper::getCommunicator());
        Dune::gridinfo(*gridPtr);
*/

        // instantiate and run the concrete problem
        Problem problem(&grid, dt, tEnd);
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
