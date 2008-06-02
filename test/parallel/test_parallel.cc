// commented lines 1454, 1464-1467 in istl/communicator.hh

#include "config.h"
#include <iostream>
#if HAVE_MPI
#include<mpi.h>
#endif
#include <dune/grid/common/gridinfo.hh>
#include <dune/common/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#include <dune/disc/operators/p1operator.hh>
#include <dune/disc/groundwater/p1groundwater.hh>

#include <dumux/io/readstarformat.cc>
#include <dumux/timedisc/timeloop.hh>
#include "parallelboxdiffusion.hh"

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
	  typedef typename Grid::Traits::LeafIndexSet IS;
	  typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,P1Layout> VM;
	  typedef typename IS::template Codim<dim>::template Partition<All_Partition>::Iterator VertexIterator;
	  
	  VM vertexMapper(grid, grid.leafIndexSet());
	  double error = 0.0;
	  
	  VertexIterator endIt = grid.leafIndexSet().template end<dim,All_Partition>();
	  VertexIterator it = grid.leafIndexSet().template begin<dim,All_Partition>();
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

    Dune::MPIHelper::instance(argc, argv);

    // define the problem dimensions  
    const int dim=2;
    typedef double NumberType; 
    if (argc != 2 && argc != 3) {
      std::cout << "usage: test_parallel dgffilename/basefilename [refinementsteps]" << std::endl;
      return 1;
    }
    int refinementSteps = 0;
    if (argc == 3) {
    	std::string arg2(argv[2]);
    	std::istringstream is2(arg2);
    	is2 >> refinementSteps;
    }
    
    // instantiate a distributed grid with overlap
    Dune::FieldVector<double,dim> length(1.0);
    Dune::FieldVector<int,dim> size(256);
    Dune::FieldVector<bool,dim> periodic(false);
    int overlap = 1;
    typedef Dune::YaspGrid<dim,dim> GridType;

    GridType grid(MPI_COMM_WORLD, length, size, periodic, overlap);

    // create a grid object
//      typedef Dune::ALUCubeGrid<dim,dim> GridType; 

// #ifdef DGF
//     // create grid pointer
//     Dune::GridPtr<GridType> gridPtr( argv[1] );
//     // grid reference 
//     GridType& grid = *gridPtr;
// #else
//     GridType grid;
//     readStarFormat(grid, argv[1]);
// #endif

//     grid.loadBalance();

    if (refinementSteps)
    	grid.globalRefine(refinementSteps);

    Dune::gridinfo(grid);
    
    Dune::Timer timer;
    timer.reset();

    DiffusionParameters<GridType,NumberType> problem;
    
    typedef Dune::LeafP1ParallelBoxDiffusion<GridType, NumberType> Diffusion;
    Diffusion diffusion(grid, problem);
    
    Dune::TimeLoop<GridType, Diffusion> timeloop(0, 1, 1, "test_parallel", 1);
    
    timeloop.execute(diffusion);

    double discreteErr = discreteError(grid, *diffusion, problem); 
    grid.comm().sum(&discreteErr, 1);
    double elapsedTime = timer.elapsed(); 
    grid.comm().max(&elapsedTime, 1); 

    if (grid.comm().rank() == 0) {
      std::cout << "discrete error = " << discreteErr << std::endl;    
      std::cout << "Calculation took " << elapsedTime << " seconds." << std::endl;
    }

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}

/*

Dune::Amg::Transfer<int, 
		    Dune::BlockVector<Dune::FieldVector<double, 1>, Dune::ISTLAllocator>, 
		    LeafP1OverlappingSchwarzCommunication<Dune::ALUCubeGrid<3, 3> > >::
prolongate(Dune::Amg::AggregatesMap<int>&, 
	   Dune::BlockVector<Dune::FieldVector<double, 1>, Dune::ISTLAllocator>&, 
	   Dune::BlockVector<Dune::FieldVector<double, 1>, Dune::ISTLAllocator>&, 
	   double, 
	   LeafP1OverlappingSchwarzCommunication<Dune::ALUCubeGrid<3, 3> >&)


Dune::Amg::Transfer<int, 
		    Dune::BlockVector<Dune::FieldVector<double, 1>, Dune::ISTLAllocator>, 
		    LeafP1OverlappingSchwarzCommunication<Dune::ALUCubeGrid<3, 3> > >::
prolongate(const Dune::Amg::AggregatesMap<V>&, 
	   Dune::BlockVector<Dune::FieldVector<double, 1>, Dune::ISTLAllocator>&, 
	   Dune::BlockVector<Dune::FieldVector<double, 1>, Dune::ISTLAllocator>&, 
	   typename V2::field_type) 

[with V1 = int, V2 = Dune::BlockVector<Dune::FieldVector<double, 1>, Dune::ISTLAllocator>, T = LeafP1OverlappingSchwarzCommunication<Dune::ALUCubeGrid<3, 3> >]
*/
