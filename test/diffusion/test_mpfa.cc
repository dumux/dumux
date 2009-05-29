// $Id$

#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>       // UGGrid
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include <dumux/material/matrixproperties.hh>
#include "dumux/diffusion/mpfa/mpfaodiffusionvelocity.hh"
#include "dumux/fractionalflow/variableclass2p.hh"

/**
 * @file
 * @brief  test case for MPFA method
 */


int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;

    // create a grid object
    typedef double NumberType; 
    typedef Dune::UGGrid<dim> GridType;
    typedef GridType::LevelGridView GridView; 

    std::stringstream dgfFileName;
    dgfFileName << "./grids/quadrilateral.dgf";

    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr( dgfFileName.str() );

    // grid reference 
    GridType& grid = *gridPtr;

    //Dune::gridinfo(grid);

    int level = 2;

    grid.globalRefine(level);

    GridView gridView(grid.levelView(level));

    //Uniform mat;
    Dune::Uniform mat;

    Dune::HomogeneousSoil<GridType, NumberType> soil(1e-10);
  
    Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, mat, mat);

    double initsat = 1;

    typedef Dune::VariableClass<GridView, NumberType> VariableType;
    VariableType variables(gridView,initsat);

    Dune::UniformProblem<GridView, NumberType, VariableType> problem(variables, mat, mat, soil, materialLaw);

    Dune::MPFAODiffusionVelocity<GridView, NumberType, VariableType> diffusion(gridView, problem);
    
    // calculate pressure on each cell
    Dune::Timer timer;
    timer.reset();
    diffusion.pressure();
    std::cout << "solving pressure took " << timer.elapsed() << " seconds" << std::endl;
    variables.vtkout("test_mpfa", level);
    
    // calculate velocity on each side
    diffusion.calcTotalVelocity();
    //printvector(std::cout, variables.velocity(), "velocity", "row", 4, 1, 3);
    
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
