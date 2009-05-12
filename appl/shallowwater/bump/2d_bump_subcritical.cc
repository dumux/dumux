//#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include "dumux/shallowwater/fvshallowwater.hh"
#include "2d_bump_subcritical_problem.hh"
#include "solidsurface_irregular.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/expliciteulerstep.hh"
#include "dumux/shallowwater/shallowvariableclass.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim = 2;

    // time loop parameters
    const double tStart = 0;
    const double tEnd = 20;
    const double cFLFactor = 0.1;
    double maxDt = 5;
    double firstDt =0.1 ;
    int modulo = 1;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(100);  N[0] = 25; N[1] = 10;
    FieldVector L(0);
    FieldVector H(0); H[0] = 25; H[1] = 10;
    GridType grid(N,L,H);

//    std::stringstream dgfFileName;
//    dgfFileName << "grids/unitcube" << GridType :: dimension << ".dgf";

    grid.globalRefine(0);

    typedef Dune::ShallowVariableClass<GridType, NumberType> VC;
    typedef Dune::SolidSurfacePlain<GridType,NumberType> Surface;

    VC variables(grid);

    Surface surface;

    typedef Dune::ShallowProblemPlain<GridType, NumberType, VC> Problem;

    Problem	problem(variables,surface,L,H);

    typedef Dune::FVShallowWater<GridType, NumberType, VC> ShallowWater;

    ShallowWater shallowWater(grid, problem);

    Dune::ExplicitEulerStep<GridType, ShallowWater> explicitEuler;

    Dune::TimeLoop<GridType, ShallowWater> timeloop(tStart, tEnd, "2dbump_subcritical", modulo, cFLFactor, maxDt, firstDt, explicitEuler);

    timeloop.execute(shallowWater);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
