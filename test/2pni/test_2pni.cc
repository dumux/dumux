//$Id$
#include "config.h"

#include "brineco2problem_2pni.hh"

#include <dune/grid/common/gridinfo.hh>

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
        typedef TTAG(BrineCO2Problem2PNI) TypeTag;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Scalar))  Scalar;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Grid))    Grid;
        typedef GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;
        typedef Dune::FieldVector<Scalar, Grid::dimensionworld> GlobalPosition;
        
        static const int dim = Grid::dimension;

        // initialize MPI, finalize is done automatically on exit
        Dune::MPIHelper::instance(argc, argv);

        // parse the command line arguments
        if (argc < 3)
            usage(argv[0]);

        // deal with the restart stuff
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

        // create the grid
        GlobalPosition lowerLeft(0.0);
        GlobalPosition upperRight;
        Dune::FieldVector<int,dim> res; // cell resolution
        upperRight[0] = 10.0;
        upperRight[1] = 10.0;
        res[0]        = 10;
        res[1]        = 10;

#if USE_UG
        ////////////////////////////////////////////////////////////
        // Make a uniform grid
        ////////////////////////////////////////////////////////////
        Grid grid;

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
                v[0] = i*(res[1]+1) + j;
                v[1] = i*(res[1]+1) + j+1;
                v[2] = (i+1)*(res[1]+1) + j;
                v[3] = (i+1)*(res[1]+1) + j+1;
                factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), v);
            }
        }

        factory.createGrid();
#else
        Grid grid(
#ifdef HAVE_MPI
                  Dune::MPIHelper::getCommunicator(),
#endif
                  upperRight, // upper right
                  res, // number of cells
                  Dune::FieldVector<bool,dim>(false), // periodic
                  2); // overlap
#endif

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
