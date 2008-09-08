// author: Jochen Fritz

#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/istl/bvector.hh>
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw_deprecated.hh"
#include "dumux/material/brookscoreylaw_deprecated.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"
#include "dumux/material/randompermeability.hh"
#include "dumux/diffusion/fv/fvdiffusion_deprecated.hh"
#include "dumux/diffusion/problems/heterogeneousproblem.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(1); N[0] = 3;
    FieldVector L(0);
    FieldVector H(300); H[0] = 600;
    GridType grid(N,L,H);

    grid.globalRefine(2);

    //Uniform mat;
    //Dune::VanGenuchtenLaw materialLaw(mat, mat);
    //Dune::BrooksCoreyLaw materialLaw(mat, mat);
    //Dune::LinearLaw materialLaw(mat, mat);

    typedef Dune::VariableClass<GridType, NumberType> VC;
    double initsat = 1;
    VC variables(grid,initsat);

    Dune::HeterogeneousProblem<GridType, NumberType, VC> problem(variables, grid, "permeab.dat", true);
    //printvector(std::cout, *(problem.permeability), "permeability", "row", 200, 1);
    //Dune::UniformProblem<GridType, NumberType> problem;
    //problem.permeability.vtkout("permeability", grid);
    //Dune::TestProblem<GridType, NumberType> problem;

    Dune::FVDiffusion<GridType, NumberType, VC> diffusion(grid, problem);
    typedef Dune::BlockVector< Dune::FieldVector<NumberType,1> > SatType;
    SatType sat(grid.levelIndexSet(0).size(0));
    int e = sat.size();
    for (int i = 0; i<e; i++) sat[i] = 0;

    diffusion.pressure();
    printvector(std::cout, problem.variables.pressure, "pressure", "row", 200, 1, 3);
    //diffusion.vtkout("pressure", 0);

    diffusion.calcTotalVelocity(0);
    printvector(std::cout, problem.variables.velocity, "velocity", "row", 4, 1, 3);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
