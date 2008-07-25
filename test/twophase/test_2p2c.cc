#include "config.h"
#include <iostream>
#ifdef HAVE_UG
#include <iomanip>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/2p2c/problems/injectionproblem.hh"
#include "dumux/2p2c/problems/uniformtwophaseproblem.hh"
#include "dumux/2p2c/fv/box2p2c.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/vangenuchtenlaw.hh"
#include "dumux/material/brookscoreylaw.hh"
#include "dumux/material/multicomponentrelations.hh"
#include "dumux/material/properties.hh"

int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions (geometry of problem)  
    const int dim=2;
    typedef double NumberType;
    Dune::FieldVector<NumberType, dim> outerLowerLeft(0);
    Dune::FieldVector<NumberType, dim> outerUpperRight(60);
    outerUpperRight[1] = 40;
    double depthBOR = 1000.0;

    // for defining e.g. a lense
    Dune::FieldVector<NumberType, dim> innerLowerLeft(4);
    innerLowerLeft[1] = 0.0;
    Dune::FieldVector<NumberType, dim> innerUpperRight(6);
    innerUpperRight[1] = 0.5;

    // count number of arguments
    if (argc != 3) {
      std::cout << "usage: test_twophase tEnd dt" << std::endl;
      return 0;
    }
    
    // define tEnd 
    std::string arg1(argv[1]);
    std::istringstream is1(arg1);
    double tEnd; is1 >> tEnd;
    // define dt
    std::string arg2(argv[2]);
    std::istringstream is2(arg2);
    double dt; is2 >> dt;


    // create a grid object
    //typedef Dune::SGrid<dim,dim> GridType; 
    //typedef Dune::YaspGrid<dim,dim> GridType; 
    typedef Dune::UGGrid<dim> GridType; 

    // use unitcube from grids (UGGrid)
    std::stringstream dgfFileName;
    dgfFileName << "dune-mux/test/twophase/grids/unitcube" 
//    dgfFileName << "grids/unitcube" 
    	<< GridType :: dimension << ".dgf";

    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr( dgfFileName.str() );

    // grid reference 
    GridType& grid = *gridPtr;

    Dune::gridinfo(grid);

    // choose fluids and properties
    Water wPhase; CO2 nPhase;
    Dune::BrooksCoreyLaw law(wPhase, nPhase);
//    Dune::LinearLaw law(wPhase, nPhase);    
    Dune::CWaterAir multicomp(wPhase, nPhase);
    //Dune::LinearLaw law(water, air);
    
    // create problem properties and geometry
    Dune::InjectionProblem<GridType, NumberType> problem(law, multicomp, 
    		outerLowerLeft, outerUpperRight, 
    		innerLowerLeft, innerUpperRight, depthBOR);

    // create two-phase two-component problem
    typedef Dune::VtkMultiWriter<GridType> MultiWriter;
    typedef Dune::Box2P2C<GridType, NumberType, MultiWriter> TwoPhaseTwoComp;
    TwoPhaseTwoComp twoPhasetwoComp(grid, problem);
    
    Dune::TimeLoop<GridType, TwoPhaseTwoComp> timeloop(0, tEnd, dt, "lens", 1);
    
    Dune::Timer timer;
    timer.reset();
    Dune::VtkMultiWriter<GridType> writer("2p2c");
    timeloop.executeMultiWriter(twoPhasetwoComp, writer);
    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;


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
  std::cout << "Please install the UG library." << std::endl;

  return 1;
}
catch (...) 
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif 
