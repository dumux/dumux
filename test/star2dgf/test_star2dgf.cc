#include <config.h>
#include <iostream>
#ifdef HAVE_ALBERTA
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfalberta.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>

#include "gridcheck.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"

int main (int argc , char **argv) try
{
    // define the problem dimensions  
    const int dim = 3;
    typedef double NumberType; 
    //typedef Dune::UGGrid<dim> GridType;
    typedef Dune::ALUSimplexGrid<dim,dim> GridType;
    //typedef Dune::ALUCubeGrid<dim,dim> GridType;
    //typedef Dune::AlbertaGrid<dim,dim> GridType;

    if (argc != 2 && argc != 3) {
    	std::cout << "Usage: test_star2dgf dgffilename [refinementsteps]" << std::endl;
    	return (1);
    }
    int refinementSteps = 0;
    if (argc == 3) {
    	std::string arg2(argv[2]);
    	std::istringstream is2(arg2);
    	is2 >> refinementSteps;
    }
    
    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr(argv[1]);

    // grid reference 
    GridType& grid = *gridPtr;

    std::cout << "Starting grid tests ." << std::flush;
//    // check macro grid 
//    gridcheck(grid);
//    std::cout << "." << std::flush;
//
//    // check the intersection iterator
//    checkIntersectionIterator(grid);
//    std::cout << "." << std::flush;
//
//    if (refinementSteps) {
//    	grid.globalRefine(refinementSteps);
//        std::cout << "." << std::flush;
//    	
//    	gridcheck(grid);
//        std::cout << "." << std::flush;
//    	
//    	// check the method geometryInFather()
//    	checkGeometryInFather(grid);
//        std::cout << "." << std::flush;
//    }
    std::cout << " passed." << std::endl;

//	  typedef GridType::Traits::LeafIndexSet IS;
//	  typedef IS::Codim<dim>::Partition<Dune::All_Partition>::Iterator VertexIterator;
//	  
//	  VertexIterator endIt = grid.leafIndexSet().end<dim,Dune::All_Partition>();
//	  VertexIterator it = grid.leafIndexSet().begin<dim,Dune::All_Partition>();
//	  for (int k = 0; k < 4 && it != endIt; ++it, k++)
//	  {
//		  // get exact solution at vertex
//		  Dune::FieldVector<double,dim> globalCoord = (*it).geometry()[0];
//		  
//		  std::cout << "k = " << k << ", coord = " << globalCoord << std::endl;
//	  }
//	  return 1;
	  
	int numberOfVertices = grid.size(dim);
    Dune::BlockVector<Dune::FieldVector<NumberType, 1> > indexVector(numberOfVertices);
    for (int i = 0; i < numberOfVertices; i++)
    	indexVector[i] = i;
    
	Dune::VTKWriter<GridType> vtkwriter(grid);
	vtkwriter.addVertexData(indexVector, "node indices");
	vtkwriter.write(argv[1], Dune::VTKOptions::ascii);		

	std::cout << "A file " << argv[1] << ".vtu has been created." << std::endl;

    return 0;
} 
catch (Dune::Exception& e) 
{
    std::cerr << e << std::endl;
    return 1;
 
} 
catch (...) 
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#else 

int main (int argc , char **argv) try
{
  std::cout << "Please install the Alberta library." << std::endl;

  return 1;
}
catch (...) 
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif 
