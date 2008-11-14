#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/istl/io.hh>
//#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include "dumux/io/vtkwriterextended.hh"
#include <dune/common/timer.hh>
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw_deprecated.hh"
#include "dumux/material/brookscoreylaw_deprecated.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"
#include "dumux/material/randompermeability.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocity_deprecated.hh"
#include "dumux/diffusion/fe/fediffusion.hh"
#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
#include "dumux/diffusion/problems/heterogeneousproblem.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
#include <dune/grid/uggrid.hh>
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

    //typedef Dune::UGGrid<dim> GridType;
    typedef Dune::SGrid<dim,dim> GridType;

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

    typedef Dune::VariableClass<GridType, NumberType> VC;
    double initsat = 1;
    VC variables(grid,initsat);

    // Dune::HeterogeneousProblem<GridType, NumberType> problem(grid, "permeab.dat", true);
    //printvector(std::cout, *(problem.permeability), "permeability", "row", 200, 1);
    Dune::UniformProblem<GridType, NumberType, VC> problem(variables);
    //Dune::AnnikaProblem<GridType, NumberType> problem;
    //problem.permeability.vtkout("permeability", grid);

    Dune::Timer timer;
    timer.reset();
    Dune::FVDiffusionVelocity<GridType, NumberType, VC> diffusion(grid, problem);
    // Dune::FEDiffusion<GridType, NumberType> diffusion(grid, problem);
    // Dune::MimeticDiffusion<GridType, NumberType, VC> diffusion(grid, problem);


    diffusion.pressure();
    std::cout << "pressure calculation took " << timer.elapsed() << " seconds" << std::endl;
    printvector(std::cout, problem.variables.pressure, "pressure", "row", 200, 1, 3);
    //diffusion.vtkout("pressure", 0);

    diffusion.calcTotalVelocity(0);
    printvector(std::cout, problem.variables.velocity, "velocity", "row", 2*dim, 1, 3);

    Dune::LeafVTKWriter<GridType> vtkwriter(grid);

    typedef Dune::BlockVector<Dune::FieldVector<double, 3> > WType;
    WType cellvelocity(grid.size(0));

    typedef Dune::BlockVector<Dune::FieldVector<double, 3> > ZType;
    ZType vertexvelocity(grid.size(dim));

     vtkwriter.faceToCellToVertex(problem.variables.velocity,cellvelocity,vertexvelocity);
     vtkwriter.addFaceData(cellvelocity,vertexvelocity,"velocity");

    vtkwriter.addCellData(problem.variables.pressure,"pressure");

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
