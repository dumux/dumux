#include <config.h>
#include <iostream>
#ifdef HAVE_UG
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/starcdreader.hh>
#include <dune/istl/bvector.hh>

#include "dumux/io/readstarformat.cc"
//#include "gridcheck.cc"
//#include "checkgeometryinfather.cc"
//#include "checkintersectionit.cc"

int main (int argc , char **argv) try
 {
     // define the problem dimensions
     const int dim = 3;
     typedef double NumberType;
     typedef Dune::UGGrid<dim> GridType;

     if (argc != 2 && argc != 3) {
         std::cout << "Usage: test_star2uggrid basefilename [refinementsteps]" << std::endl;
         return (1);
     }
     int refinementSteps = 0;
     if (argc == 3) {
         std::string arg2(argv[2]);
         std::istringstream is2(arg2);
         is2 >> refinementSteps;
     }

     GridType* grid;

     //readStarFormat(grid, argv[1]);
     Dune::StarCDReader<GridType> reader;
     grid = reader.read(argv[1]);

     std::cout << "Starting grid tests ." << std::flush;
     // check macro grid
     //gridcheck(grid);
     std::cout << "." << std::flush;

     // check the intersection iterator
     //checkIntersectionIterator(grid);
     std::cout << "." << std::flush;

     if (refinementSteps) {
         grid->globalRefine(refinementSteps);
         std::cout << "." << std::flush;

         //gridcheck(grid);
         std::cout << "." << std::flush;

         // check the method geometryInFather()
         //checkGeometryInFather(grid);
         std::cout << "." << std::flush;
     }
     std::cout << " passed." << std::endl;

     int numberOfVertices = grid->size(dim);
     Dune::BlockVector<Dune::FieldVector<NumberType, 1> > indexVector(numberOfVertices);
     for (int i = 0; i < numberOfVertices; i++)
         indexVector[i] = i;

     Dune::VTKWriter<GridType::LeafGridView> vtkwriter(grid->leafView());
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
     std::cout << "Please install the UG library." << std::endl;

     return 1;
 }
 catch (...)
 {
     std::cerr << "Generic exception!" << std::endl;
     return 2;
 }
#endif
