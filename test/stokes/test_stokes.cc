#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include <dune/disc/stokes/dgstokes.hh>
#include <dune/disc/stokes/l2error.hh>

int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;
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
    Example<dim, NumberType> exactSolution;
    DirichletBoundary<GridType> dirichletBoundary(exactSolution);
    RightHandSide<GridType> rightHandSide(exactSolution);
    
	typedef Dune::DGStokes<GridType, vOrder, pOrder> DGStokes;
	DGStokes dGStokes(grid, exactSolution, parameters, dirichletBoundary, rightHandSide, refinementSteps); 
	dGStokes.assembleStokesSystem();
	dGStokes.solveStokesSystem();
	
	std::cout << "L2Error: ";
	for (int i = 0; i < dim+1; i++)
		std::cout << dGStokes.l2errorStokesSystem(i) << ", ";
	std::cout << std::endl;
	
	dGStokes.vtkout("test_stokes", 0);

	return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
