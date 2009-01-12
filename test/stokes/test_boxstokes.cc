#include "config.h"
#include <iostream>
#define DUMMY
#ifdef DUMMY
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/timedisc/timeloop.hh"
#include "boxstokes.hh"
#include "yxproblem.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim = 2;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;

    if (argc != 2 && argc != 3) {
        std::cout << "Usage: test_boxstokes dgffilename [refinementsteps]" << std::endl;
        return (1);
    }
    int refinementSteps = 0;
    if (argc == 3) {
        std::string arg2(argv[2]);
        std::istringstream is2(arg2);
        is2 >> refinementSteps;
    }

    Dune::GridPtr<GridType> gridPtr( argv[1] );
    GridType& grid = *gridPtr;

    if (refinementSteps)
        grid.globalRefine(refinementSteps);

    Dune::YXProblem<GridType, double> problem;
    typedef Dune::LeafP1BoxStokes<GridType, NumberType, dim> BoxStokes;
    BoxStokes boxStokes(grid, problem);

    Dune::TimeLoop<GridType, BoxStokes> timeloop(0, 1, 1, "test_boxstokes", 1);

    timeloop.execute(boxStokes);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
#else

int main (int argc , char **argv) try
{
  std::cout << "This test is not finished yet." << std::endl;

  return 1;
}
catch (...)
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif
