#include "config.h"

#include "new_blobproblem.hh"

#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/uggrid.hh>

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
        const int dim = 2;
        typedef double                                Scalar;
        typedef Dune::UGGrid<dim>                     Grid;
        typedef Dune::NewBlobProblem<Grid, Scalar>    Problem;
        typedef Problem::DomainTraits::GlobalPosition GlobalPosition;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments
        if (argc < 3)
            usage(argv[0]);

        int argPos = 1;
        bool restart = false;
        double restartTime = 0;
        if (std::string("--restart") == argv[argPos]) {
            restart = true;
            ++argPos;

            std::istringstream(argv[argPos++]) >> restartTime;
        }

        if (argc - argPos != 2) {
            usage(argv[0]);
        }

        double tEnd, dt;
        std::istringstream(argv[argPos++]) >> tEnd;
        std::istringstream(argv[argPos++]) >> dt;

        ////////////////////////////////////////////////////////////
        // Make a uniform grid
        ////////////////////////////////////////////////////////////
        Grid grid;

        GlobalPosition upperRight;
        Dune::FieldVector<int,dim> res; // cell resolution
        upperRight[0] = 300.0;
        res[0]        = 20;
        upperRight[1] = 300.0;
        res[1]        = 20;

        Dune::GridFactory<Grid> factory(&grid);
        for (int i=0; i<=res[0]; i++) {
            for (int j=0; j<=res[1]; j++) {
                Dune::FieldVector<double,2> pos;
                pos[0] = upperRight[0]*double(i)/res[0];
                pos[1] = upperRight[1]*double(j)/res[1];
                factory.insertVertex(pos);
            }
        }

        for (int i=0; i<res[0]; i++) {
            for (int j=0; j<res[1]; j++) {
                std::vector<unsigned int> v(4);
                v[0] = i*(res[0]+1) + j;
                v[1] = i*(res[0]+1) + j+1;
                v[2] = (i+1)*(res[0]+1) + j;
                v[3] = (i+1)*(res[0]+1) + j+1;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), v);
            }
        }

        factory.createGrid();

        ////////////////////////////////////////////////////////////
        // instantiate and run the concrete problem
        ////////////////////////////////////////////////////////////
        Problem problem(&grid, dt, tEnd);

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
        throw;
    }

    return 3;
}
