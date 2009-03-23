#include "config.h"
#include "nimasdiffusionproblem.hh"

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
//#include <dune/grid/io/file/dgfparser/dgfug.hh>

#include <dune/common/exceptions.hh>
#include <dune/common/mpihelper.hh>

#include <iostream>
#include <boost/format.hpp>

int main(int argc, char** argv)
{
    try {
        const int dim = 2;
        // Set the type for scalar values (should be one of float, double
        // or long double)
        typedef double                            Scalar;
        typedef Dune::SGrid<dim,dim>              Grid;
        typedef Dune::NimasDiffusionProblem<Grid, Scalar>  Problem;

        //    For defining a SGrid ////////////////////////////////////////
        typedef Dune::SGrid<dim,dim> Grid;
        typedef Dune::FieldVector<Grid::ctype,dim> FieldVector;
        Dune::FieldVector<int,dim> N(4); N[1]=25; //number of cells
        FieldVector L(0); // lower left corner
        FieldVector H(0.073); H[1]=0.25; // upper right corner
        Grid grid(N,L,H);
        ////////////
//        typedef Dune::GridPtr<Grid>               GridPointer;


        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments for the program
        if (argc != 3) {
//            std::cout << boost::format("usage: %s gridFile.dgf tEnd dt\n")%argv[0];
            std::cout << boost::format("usage: %s tEnd dt\n")%argv[0];
            return 1;
        }
        double tEnd, dt;
//        const char *dgfFileName = argv[1];
        std::istringstream(argv[1]) >> tEnd;
        std::istringstream(argv[2]) >> dt;

        Dune::gridinfo(grid);

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
