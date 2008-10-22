#include "config.h"
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>



int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;

    // create a grid object
    typedef Dune::UGGrid<dim> GridType; 
    typedef GridType::ctype DT;
    typedef double NumberType; 

    std::stringstream dgfFileName;
    dgfFileName << "grids/mesh2_1.dgf";

    // create grid pointer
    Dune::GridPtr<GridType> gridPtr1( dgfFileName.str() );    
    Dune::GridPtr<GridType> gridPtr2( dgfFileName.str() );

    // grid reference
    GridType& grid1 = *gridPtr1;
    GridType& grid2 = *gridPtr2;

    Dune::gridinfo(grid1);
    Dune::gridinfo(grid2);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
