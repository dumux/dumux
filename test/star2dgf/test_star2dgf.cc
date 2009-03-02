#include <config.h>
#include <iostream>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/bvector.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>

int main (int argc , char **argv) try
 {
     const int dim = 3;
     typedef double NumberType;
     typedef Dune::UGGrid<dim> GridType;
     //typedef Dune::ALUCubeGrid<dim,dim> GridType;

     if (argc != 2) {
         std::cout << "Usage: test_dgf dgffilename" << std::endl;
         return (1);
     }

     Dune::GridPtr<GridType> gridPtr(argv[1]);

     GridType& grid = *gridPtr;

     int numberOfVertices = grid.size(dim);
     Dune::BlockVector<Dune::FieldVector<NumberType, 1> > indexVector(numberOfVertices);
     for (int i = 0; i < numberOfVertices; i++)
         indexVector[i] = i;

     Dune::VTKWriter<GridType::LeafGridView> vtkwriter(grid.leafView());
     vtkwriter.addVertexData(indexVector, "node indices");
     vtkwriter.write(argv[1], Dune::VTKOptions::ascii);

     //     std::cout << "A file " << argv[1] << ".vtu has been created." << std::endl;

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
