#include "config.h"
#include <iostream>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include "dumux/operators/p1operatorextended.hh"
#include <dune/disc/groundwater/p1groundwater.hh>

#include <dumux/io/readstarformat.cc>
#include <dumux/timedisc/timeloop.hh>
#include "boxdiffusion.hh"

#define DGF

namespace Dune
{
template<int dim>
struct P1Layout
{
	bool contains (Dune::GeometryType gt)
	{
		return gt.dim() == 0;
	}
}; 
template<class Grid, class Solution, class Problem> 
double discreteError(const Grid& grid, const Solution& solution, const Problem& problem)
{
	  enum{dim=Grid::dimension};
		typedef typename Grid::LeafGridView GV;
	    typedef typename GV::IndexSet IS;
	  typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,P1Layout> VM;
		typedef typename GV::template Codim<dim>::Iterator VertexIterator;
	  
	  VM vertexMapper(grid, grid.leafIndexSet());
	  double error = 0.0;
	  const GV& gridview(grid.leafView());
	  
	  VertexIterator endIt = gridview.template end<dim>();
	  VertexIterator it = gridview.template begin<dim>();
	  for (; it != endIt; ++it)
	  {
		  // get exact solution at vertex
		  FieldVector<double,dim> globalCoord = (*it).geometry()[0];
		  double exact = problem.exact(globalCoord);

		  // get approximate solution at vertex
		  int globalId = vertexMapper.map(*it);
		  double approximate = (*solution)[globalId];
		  
		  error += (exact - approximate)*(exact - approximate);
	  }
		  
	  return sqrt(error);
}
}

int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;
    typedef double NumberType; 
    if (argc != 2 && argc != 3) {
      std::cout << "usage: box3d dgffilename/basefilename [refinementsteps]" << std::endl;
      return 1;
    }
    int refinementSteps = 0;
    if (argc == 3) {
    	std::string arg2(argv[2]);
    	std::istringstream is2(arg2);
    	is2 >> refinementSteps;
    }
    
    // create a grid object
    typedef Dune::SGrid<dim,dim> GridType; 
    //typedef Dune::YaspGrid<dim,dim> GridType; 
    //typedef Dune::UGGrid<dim> GridType; 
    //typedef Dune::ALUCubeGrid<dim,dim> GridType; 

#ifdef DGF
    // create grid pointer
    Dune::GridPtr<GridType> gridPtr( argv[1] );
    // grid reference 
    GridType& grid = *gridPtr;
#else
    GridType grid;
    readStarFormat(grid, argv[1]);
#endif
    
    if (refinementSteps)
    	grid.globalRefine(refinementSteps);

    Dune::gridinfo(grid);
    
    DiffusionParameters<GridType,NumberType> problem;
    
    typedef Dune::LeafP1BoxDiffusion<GridType, NumberType> Diffusion;
    Diffusion diffusion(grid, problem);
    
    Dune::TimeLoop<GridType, Diffusion> timeloop(0, 1, 1, "box3d", 1);
    
    timeloop.execute(diffusion);

	std::cout << "discrete error = " << discreteError(grid, *diffusion, problem) << std::endl;
	return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
