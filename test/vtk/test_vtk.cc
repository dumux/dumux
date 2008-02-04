#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include "dumux/io/vtkwriterextended.hh"
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
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
 
int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;

    // create a grid object
    typedef double NumberType; 
    typedef Dune::SGrid<dim,dim> GridType; 
    //typedef Dune::YaspGrid<dim,dim> GridType; 
    //typedef Dune::UGGrid<dim> GridType; 

    // use unitcube from grids 
    std::stringstream dgfFileName;
    dgfFileName << "grids/unitcube" << GridType :: dimension << ".dgf";

    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr( dgfFileName.str() );

    // grid reference 
    GridType& grid = *gridPtr;

    Dune::gridinfo(grid);

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
    Dune::FVDiffusion<GridType, NumberType> diffusion(grid, problem);
    //Dune::FEDiffusion<GridType, NumberType> diffusion(grid, problem);
    //Dune::MimeticDiffusion<GridType, NumberType> diffusion(grid, problem);
    
    
    diffusion.pressure();
    std::cout << "pressure calculation took " << timer.elapsed() << " seconds" << std::endl;
    printvector(std::cout, *diffusion, "pressure", "row", 200, 1, 3);
    //diffusion.vtkout("pressure", 0);
    
    const int blocksize = 2*dim;
    typedef Dune::FieldVector<double, dim> R1;
    typedef Dune::BlockVector< Dune::FieldVector<R1, blocksize> > VType;
    VType velocity(grid.size(0));
    diffusion.totalVelocity(velocity);
    printvector(std::cout, velocity, "velocity", "row", 4, 1, 3);
    
	//Dune::VTKWriterExtended<GridType, GridType::Codim<0>::LevelIndexSet> 
	//	vtkwriter(grid, grid.levelIndexSet( grid.maxLevel() ));
	Dune::VTKWriter<GridType> vtkwriter(grid);
	//vtkwriter.addFaceData(velocity, "face velocities");
	vtkwriter.addCellData(*diffusion, "pressure");
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
