#include "config.h"
#include <iostream>
#ifdef HAVE_UG
#include <iomanip>
//#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
//#include <dune/grid/io/file/dgfparser/dgfalu.hh>
//#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
//#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
//#include "dumux/twophase/problems/uniformtwophaseproblem.hh"
#include "waterco2problem.hh"
//#include "dumux/material/properties.hh"
#include "dumux/2p2cni/fv/boxco2.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/io/readstarformat.cc"
int main(int argc, char** argv) 
{
  try{
    // define the problem dimensions  
    const int dim=2;
    typedef double NumberType; 
    Dune::FieldVector<NumberType, dim> outerLowerLeft(0);
    Dune::FieldVector<NumberType, dim> outerUpperRight(4);
    Dune::FieldVector<NumberType, dim> innerLowerLeft(1);
    innerLowerLeft[1] = 2;
    Dune::FieldVector<NumberType, dim> innerUpperRight(4);
    innerUpperRight[1] = 3;
    if (argc != 4) {
      std::cout << "usage: 2p2cni basefilename tEnd dt" << std::endl;
      return 0;
    }
    	std::string arg1(argv[2]);
	std::istringstream is1(arg1);
	double tEnd;
	is1 >> tEnd;
	std::string arg2(argv[3]);
	std::istringstream is2(arg2);
	double dt;
	is2 >> dt;

    // create a grid object
    typedef Dune::SGrid<dim,dim> GridType; 
    //typedef Dune::YaspGrid<dim,dim> GridType; 
    //typedef Dune::UGGrid<dim> GridType; 
	//typedef Dune::ALUSimplexGrid<dim,dim> GridType; 
	
    Dune::GridPtr<GridType> gridPointer(argv[1]);
    GridType& grid = *gridPointer;
    //readStarFormat(grid, argv[1]);
    //grid.createLGMGrid(argv[1]);

     Dune::gridinfo(grid);
    
    Brine wPhase(0.2);
    CO2 nPhase(0.05);

//     Dune::LinearLaw law(brine, co2);
//     Dune::CO2Problem2D<GridType, NumberType> problem(law, 1.e7); 
      Dune::BrooksCoreyLaw law(wPhase, nPhase);
      Dune::CBrineCO2 multicomp(wPhase, nPhase);

      Dune::WaterCO2Problem<GridType, NumberType> problem(law, multicomp, 9.548355e6, 310, 0.2, 0.05, 10000, 2); 

      typedef Stupid::VtkMultiWriter<GridType> MultiWriter;
    typedef Dune::BoxCO2<GridType, NumberType, MultiWriter> TwoPhase;
    TwoPhase twoPhase(grid, problem);


    
//    Dune::TimeLoop<GridType, TwoPhase> timeloop(0, tEnd, dt, 1.e3, "co2", 5);
    Dune::TimeLoop<GridType, TwoPhase> timeloop(0, tEnd, dt, "co2", 5);
    
    Dune::Timer timer;
    timer.reset();
    timeloop.execute(twoPhase);
    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;
     
	//std::cout << twoPhase.injected() << " kg CO2 injected." << std::endl;

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
