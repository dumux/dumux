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
//    try {

        // Set the type for scalar values (should be one of float, double
        // or long double)
        const int dim = 2;
        typedef double                                Scalar;
//        typedef Dune::ALUSimplexGrid<dim,dim>         Grid;
#define USE_UG 0
#if USE_UG
        typedef Dune::UGGrid<dim>                     Grid;
#else
        typedef Dune::YaspGrid<dim>                   Grid;
#endif
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

#if USE_UG
        // //////////////////////////////////////////////////////////
        //   Make a uniform grid on rank 0 for testing
        // //////////////////////////////////////////////////////////
        int n = 11;
        Grid grid;
        Dune::GridFactory<Grid> factory(&grid);

        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                Dune::FieldVector<double,2> pos;
                pos[0] = 300*double(i)/(n-1);
                pos[1] = 300*double(j)/(n-1);
                factory.insertVertex(pos);
            }
        }
    
        for (int i=0; i<n-1; i++) {
            for (int j=0; j<n-1; j++) {
                std::vector<unsigned int> v(4);
                v[0] = i*n + j;
                v[1] = i*n + j+1;
                v[2] = (i+1)*n + j;
                v[3] = (i+1)*n + j+1;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), v);
            }
        }

        //std::auto_ptr<Dune::UGGrid<2> > grid(factory.createGrid());
        factory.createGrid();
#else

        GlobalPosition upperRight;
        Dune::FieldVector<int,dim> res; // cell resolution
        upperRight[0] = 300.0;
        res[0]        = 6;

        upperRight[1] = 300.0;
        res[1]        = 6;

        Grid grid(
#ifdef HAVE_MPI
                  Dune::MPIHelper::getCommunicator(),
#endif
                  upperRight, // upper right
                  res, // number of cells
                  Dune::FieldVector<bool,dim>(false), // periodic
                  2); // overlap
#endif
        
        // instantiate and run the concrete problem
        Problem problem(&grid, dt, tEnd);
        if (!problem.simulate())
            return 2;
        
        return 0;
/*    }
    catch (Dune::Exception &e) {
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!\n";
        throw;
    }
*/
    return 3;
}
