#include "config.h"
#include <iostream>
#include <iomanip>
//#ifdef HAVE_UG
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/sgrid.hh>
//#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw.hh"
#include "dumux/material/brookscoreylaw.hh"
#include "dumux/material/vangenuchtenlaw.hh"
#include "dumux/material/randompermeability.hh"
#include "dumux/diffusion/fv/fvdiffusion.hh"
//#include "dumux/diffusion/fe/fediffusion.hh"
//#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
#include "dumux/diffusion/problems/heterogeneousproblem.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
 
int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;

    // create a grid object
    typedef double NumberType; 
    typedef Dune::SGrid<dim,dim> GridType; 
    //typedef Dune::ALUSimplexGrid<dim,dim> GridType; 
    //typedef Dune::YaspGrid<dim,dim> GridType; 
    //typedef Dune::UGGrid<dim> GridType; 

    // use unitcube from grids 
    std::stringstream dgfFileName;
    //dgfFileName << "grids/skew.dgf";
    dgfFileName << "grids/unitcube" << GridType :: dimension << ".dgf";

    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr( argv[1] );

    // grid reference 
    //GridType& grid = *gridPtr;
    Dune::FieldVector<GridType::ctype,dim> L(0);
    Dune::FieldVector<GridType::ctype,dim> R(300);
    Dune::FieldVector<int,dim> N(10);           
    GridType grid(N,L,R);

    grid.globalRefine(3);

    Dune::gridinfo(grid);

    Dune::SimpleProblem<GridType, NumberType> satprob;
    //Uniform mat;
    //Dune::VanGenuchtenLaw materialLaw(mat, mat);
    //Dune::BrooksCoreyLaw materialLaw(mat, mat);
    //Dune::LinearLaw materialLaw(mat, mat);
    
    //Dune::HeterogeneousProblem<GridType, NumberType> problem(grid, "permeab.dat", true);
    //printvector(std::cout, *(problem.permeability), "permeability", "row", 200, 1);
    Dune::UniformProblem<GridType, NumberType> problem;
    //problem.permeability.vtkout("permeability", grid);

    Dune::Timer timer;
    timer.reset();
    //Dune::FEDiffusion<GridType, NumberType> diffusion(grid, problem);
    Dune::FVDiffusion<GridType, NumberType> diffusion(grid, problem, satprob, grid.maxLevel());
    //Dune::MimeticDiffusion<GridType, NumberType> diffusion(grid, problem, satprob, grid.maxLevel());
    
    
    diffusion.pressure();
    std::cout << "pressure calculation took " << timer.elapsed() << " seconds" << std::endl;
    printvector(std::cout, *diffusion, "pressure", "row", 200, 1, 3);
    diffusion.vtkout("fv", 0);
    
    //const int blocksize = 2*dim;
    //typedef Dune::FieldVector<double, dim> R1;
    //typedef Dune::BlockVector< Dune::FieldVector<R1, blocksize> > VType;
    //VType velocity(grid.size(0));
    //diffusion.totalVelocity(velocity);
    //printvector(std::cout, velocity, "velocity", "row", 4, 1, 3);
    
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
//#else 
//
//int main (int argc , char **argv) try
//{
//  std::cout << "Please install the UG library." << std::endl;
//
//  return 1;
//}
//catch (...) 
//{
//    std::cerr << "Generic exception!" << std::endl;
//    return 2;
//}
//#endif 
