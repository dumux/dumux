
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <iomanip>
#include"dune/common/mpihelper.hh" // An initializer of MPI
#include"dune/common/exceptions.hh" // We use exceptions
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
//#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/dgfparser.hh>
//#include <dune/grid/io/file/dgfparser/dgfug.hh>
//#include <dune/grid/io/file/dgfparser/dgfs.hh>
//#include <dune/grid/io/file/dgfparser/dgfalu.hh>
//#include <dune/grid/io/file/dgfparser/dgfalberta.hh>
//#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
//#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include <boost/format.hpp>
//#include "richardsproblemPw_couple.hh"
//#include "boxpw_couple.hh"
//#include "./timedisc/timeloop.hh"
//#include "dumux/material/phaseproperties/phaseproperties2p.hh"
//#include "dumux/material/twophaserelations.hh"
//#include "dune_multiscale_headers.hh"
#include "new_richardsproblempipe.hh"

int main(int argc, char** argv)
{
    try{
        // Initialize MPI
        Dune::MPIHelper::instance(argc, argv);

        // set up the grid
        const int dim = 3;
        typedef double                             Scalar;
        typedef Dune::ALUCubeGrid<dim,dim>         Grid;
        //        typedef Dune::ALUSimplexGrid<dim,dim>     Grid;
        //        typedef Dune::UGGrid<dim>                 Grid;
        typedef Dune::GridPtr<Grid>          GridPointer;

        if (argc != 5) {
            std::cout << "usage: pipe_richards_incomp_couple_timeloop dgffilename/basefilename tEnd dt dtmax" << std::endl;
            return 1;
        }
        double tEnd, dt, dtmax;
        const char *dgfFileName = argv[1];
        std::istringstream(argv[2]) >> tEnd;
        std::istringstream(argv[3]) >> dt;
        std::istringstream(argv[4]) >> dtmax;

        // create grid pointer, Grid is defined by gridtype.hh
        GridPointer gridPtr( dgfFileName );

        // grid reference
        Grid& grid = *gridPtr;
        Dune::gridinfo(grid);

        typedef Dune::NewRichardsProblemPipe<Grid, Scalar>   Problem;

        // instantiate and run the concrete problem
        Problem problem(gridPtr, dt, dtmax, tEnd);
        if (!problem.simulate())
            return 2;

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}

