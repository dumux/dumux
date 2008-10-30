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
#include "dumux/stokes/dgstokes.hh"
#include "dumux/stokes/l2error.hh"
#include "dumux/stokes/h1error.hh"
#include "yxproblem.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim = 2;
    const int vOrder = 2;
    const int pOrder = 1;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;

    if (argc != 2 && argc != 3) {
    	std::cout << "Usage: test_stokes dgffilename [refinementsteps]" << std::endl;
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

    DGStokesParameters parameters;
    Dune::YXProblem<GridType, double> problem;
    typedef Dune::DGStokes<GridType, vOrder, pOrder> DGStokes;
    DGStokes dGStokes(grid, problem, parameters);
    dGStokes.assembleStokesSystem();

//    printmatrix(std::cout, dGStokes.matrix(), "stiffness matrix", "row", 11, 4);
//    printvector(std::cout, dGStokes.rhs(), "right hand side", "row", 200, 1, 3);
    dGStokes.solveStokesSystem();
    dGStokes.vtkout("test_stokes", 0);
    printvector(std::cout, dGStokes.sol(), "solution", "row", 200, 1, 3);

    //	std::cout << "L2Error velocity: ";
    //	for (int i = 0; i < dim; i++)
    //		std::cout << dGStokes.l2errorStokesSystem(i) << ", ";
    //	std::cout << std::endl;
    //	std::cout << "L2Error pressure: "<< dGStokes.l2errorStokesSystem(dim) << std::endl;
    //	std::cout << "H1Error velocity: ";
    //	for (int i = 0; i < dim; i++)
    //		std::cout << dGStokes.h1errorStokesSystem(i) << ", ";
    //	std::cout << std::endl;

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
