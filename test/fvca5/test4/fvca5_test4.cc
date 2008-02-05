#include "config.h"
#include <iostream>
#ifdef HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/fe/fediffusion.hh"
#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
#include "fvca5test4problem.hh"
#include "../benchmarkresult.hh"

int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;

    // create a grid object
    typedef double NumberType; 
    typedef Dune::UGGrid<dim> GridType; 

    if (argc != 2 && argc != 3) {
    	std::cout << "Usage: fvca5_test4 dgffilename [refinementsteps]" << std::endl;
    	return (1);
    }
    int refinementSteps = 0;
    if (argc == 3) {
    	std::string arg2(argv[2]);
    	std::istringstream is2(arg2);
    	is2 >> refinementSteps;
    }
    
    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr( argv[1] );

    // grid reference 
    GridType& grid = *gridPtr;

    if (refinementSteps)
    	grid.globalRefine(refinementSteps);

	Dune::FVCA5Test4Problem<GridType, NumberType> problem;
    Dune::SimpleProblem<GridType, NumberType> satprob;

    Dune::Timer timer;
    timer.reset();
    //Dune::FEDiffusion<GridType, NumberType> diffusion(grid, problem);
    //Dune::FVDiffusion<GridType, NumberType> diffusion(grid, problem);
    Dune::MimeticDiffusion<GridType, NumberType> diffusion(grid, problem, satprob, grid.maxLevel());
    
    diffusion.pressure();
    std::cout << "pressure calculation took " << timer.elapsed() << " seconds" << std::endl;
    //printvector(std::cout, *diffusion, "pressure", "row", 200, 1, 3);
    
    Dune::BenchmarkResult result;
    result.evaluate(grid, problem, diffusion);  
    std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
    std::cout.setf(std::ios_base::uppercase);
    std::cout.precision(2);
    
    std::cout << "sumflux = flux0 + flux1 + fluy0 + fluy1 - sumf \n        = " 
    	<< result.flux0 << " + " << result.flux1 << " + " 
    	<< result.fluy0 << " + " << result.fluy1 << " - " 
    	<< result.sumf << "\n        = " << result.sumflux << std::endl;
    std::cout << "energy error eren = |ener1 - ener2|/max(ener1,ener2) = |"
    	<< result.ener1 << " - " << result.ener2 << "|/" << std::max(result.ener1, result.ener2)
    	<< " = " << result.eren << std::endl;
    std::cout << "umin = " << result.uMin << std::endl;
    std::cout << "umax = " << result.uMax << std::endl;
    
    
    diffusion.vtkout("fvca5_test4", 0);


    
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
  std::cout << "Please install the UG library." << std::endl;

  return 1;
}
catch (...) 
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif 
