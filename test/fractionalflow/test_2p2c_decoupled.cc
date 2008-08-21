#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/phaseproperties.hh"
#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>
#include "dumux/transport/problems/testproblem_2p2c.hh"
#include "timeloop.hh"
#include "decoupled2p2c.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include "dumux/timedisc/expliciteulerstep.hh"
 
int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType; 
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector; 
    Dune::FieldVector<int,dim> N(2); N[0] = 2;                   
    FieldVector L(0); 
    FieldVector H(300); H[0] = 300; 
    GridType grid(N,L,H);  
 
    grid.globalRefine(0);
    
    double tStart = 0; 
    double tEnd = 1e5;
    char* fileName("timeloop");
    int modulo = 1; 
    double cFLFactor = 0.99;
    
    Dune::L_air_water wetmat;
    Dune::G_air_water nonwetmat;
    Dune::Homogeneoussoil<GridType, NumberType> soil;
    Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wetmat, nonwetmat);
    
    Dune::VariableClass2p2c<GridType,NumberType> var(grid);
    
    typedef Dune::Testproblem_2p2c<GridType, NumberType> TransProb;
    TransProb problem(grid, var, wetmat, nonwetmat, soil, grid.maxLevel(), materialLaw, false);
//    problem.permeability.vtkout("permeability", grid);

    typedef Dune::Decoupled2p2c<GridType, NumberType> ModelType;
    ModelType model(grid, problem, grid.maxLevel());

    Dune::ExplicitEulerStep<GridType, ModelType> timestep;
    Dune::TimeLoop<GridType, ModelType > timeloop(tStart, tEnd, fileName, modulo, cFLFactor, 1e100, 1e100, timestep);
    
    Dune::Timer timer;
    timer.reset();
    timeloop.execute(model);
    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;
    
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
