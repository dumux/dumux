#include "config.h"
#include <iostream>
#define DUMMY
#ifdef DUMMY
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/disc/functions/p1function.hh>

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;
    typedef double NumberType;

//    Dune::FieldVector<double,dim> length(1.0);
//    length[0] = 2.0;
//    Dune::FieldVector<int,dim> size(1);
//    size[0] = 2*size[1];
//    Dune::FieldVector<bool,dim> periodic(false);
//    int overlap = 0;
//    typedef Dune::YaspGrid<dim,dim> GridType;
//    GridType grid(length,size,periodic,overlap);

    typedef Dune::UGGrid<dim> GridType;
    //typedef Dune::SGrid<dim,dim> GridType;
    Dune::GridPtr<GridType> gridPtr(argv[1]);
    GridType& grid(*gridPtr);

    typedef Dune::LeafP1Function<GridType, NumberType> P1Function;
    P1Function u(grid);
    (*u)[0] = 0;
    (*u)[1] = 0;
    (*u)[2] = 0;
    (*u)[3] = 0;
    (*u)[4] = 1.0;
    (*u)[5] = 2.0;

    Dune::FieldVector<double, dim> global(0);
    std::cout << global << ": " << u.eval(0, global) << std::endl;
    global[0] = 0.5; global[1] = 0.5;
    std::cout << global << ": " << u.eval(0, global) << std::endl;
    global[0] = 1.5; global[1] = 0.5;
    std::cout << global << ": " << u.eval(0, global) << std::endl;
    global[0] = 0.33; global[1] = 0.67;
    std::cout << global << ": " << u.eval(0, global) << std::endl;

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
