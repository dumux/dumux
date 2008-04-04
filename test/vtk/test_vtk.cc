#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/io/vtkwriterextended.hh"
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw.hh"
#include "dumux/material/brookscoreylaw.hh"
#include "dumux/material/vangenuchtenlaw.hh"
#include "dumux/material/randompermeability.hh"
#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/fe/fediffusion.hh"
#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
#include "dumux/diffusion/problems/heterogeneousproblem.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
//#include "annikaproblem.hh" 
 
int main(int argc, char** argv) 
{
  if(argc<2)
      std::cout << "missing argument: dgfFileName" << std::endl;

  try{
    // define the problem dimensions  
    const int dim=2;

    // create a grid object
    typedef double NumberType; 
  
    typedef Dune::UGGrid<dim> GridType; 
    //typedef Dune::ALUSimplexGrid<dim,dim> GridType;
    
    // use unitcube from grids 
    std::stringstream dgfFileName;

    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr(argv[1]);
      
    // grid reference 
    GridType& grid = *gridPtr;

    Dune::gridinfo(grid);

    //Uniform mat;
    //Dune::VanGenuchtenLaw materialLaw(mat, mat);
    //Dune::BrooksCoreyLaw materialLaw(mat, mat);
    //Dune::LinearLaw materialLaw(mat, mat);
    
    // Dune::HeterogeneousProblem<GridType, NumberType> problem(grid, "permeab.dat", true);
    //printvector(std::cout, *(problem.permeability), "permeability", "row", 200, 1);
    Dune::UniformProblem<GridType, NumberType> problem;
    //Dune::AnnikaProblem<GridType, NumberType> problem;
    //problem.permeability.vtkout("permeability", grid);

    Dune::Timer timer;
    timer.reset();
    //  Dune::FVDiffusion<GridType, NumberType> diffusion(grid, problem);
    // Dune::FEDiffusion<GridType, NumberType> diffusion(grid, problem);
     Dune::MimeticDiffusion<GridType, NumberType> diffusion(grid, problem);
    
    
    diffusion.pressure();
    std::cout << "pressure calculation took " << timer.elapsed() << " seconds" << std::endl;
    printvector(std::cout, *diffusion, "pressure", "row", 200, 1, 3);
    //diffusion.vtkout("pressure", 0);
    
    const int blocksize = 2*dim;
    typedef Dune::FieldVector<NumberType, dim> R1;
    typedef Dune::BlockVector< Dune::FieldVector<R1, blocksize> > VType;
    VType velocity(grid.size(0));
    diffusion.totalVelocity(velocity);
    printvector(std::cout, velocity, "velocity", "row", 2*dim, 1, 3);

    Dune::VTKWriter<GridType> vtkwriter(grid);
   
    typedef Dune::BlockVector<R1> WType;
    WType cellvelocity(grid.size(0));

    typedef Dune::BlockVector<R1> ZType;
    ZType vertexvelocity(grid.size(dim));

    vtkwriter.faceToCellToVertex(velocity,cellvelocity,vertexvelocity);
    vtkwriter.addFaceData(cellvelocity,vertexvelocity,"velocity");

    vtkwriter.addCellData(*diffusion,"pressure");
  
    vtkwriter.write("test_vtk", Dune::VTKOptions::ascii);		

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
