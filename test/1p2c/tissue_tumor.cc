//$Id$

#include "config.h"
#include <iostream>
#define DUMMY
//#undef DUMMY
#ifdef DUMMY
//#ifdef HAVE_UG
#include <iomanip>
//#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
//#include <dune/grid/io/file/dgfparser/dgfug.hh>
//#include <dune/grid/io/file/dgfparser/dgfalu.hh>
//#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
//#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "tissue_tumor_problem.hh"
#include "tissue_soilproperties.hh"

#include "dumux/material/phaseproperties/phaseproperties1p.hh"
#include "dumux/material/property_baseclasses.hh"

#include "dumux/1p2c/fv/box1p2c.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/io/readstarformat.cc"
//#include "dumux/io/vtkmultiwriter.hh"


int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;
    typedef double NumberType;


    if (argc != 4) {
      std::cout << "usage: 1p2c basefilename tEnd dt" << std::endl;
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

     // choose fluids and properties
     Dune::InterstitialFluid phase;
     Dune::TissueSoil<GridType, NumberType> soil;


     Dune::TissueTumorProblem<GridType, NumberType> problem(phase, soil);

 

     typedef Dune::VtkMultiWriter<GridType::LeafGridView> MultiWriter;
     typedef Dune::Box1P2C<GridType, NumberType, MultiWriter> OnePhaseTwoComp;
     OnePhaseTwoComp onePhasetwoComp(grid, problem);

     Dune::TimeLoop<GridType, OnePhaseTwoComp, false> timeloop(0, tEnd, dt, "tissue_tumor", 1);

//      Dune::Timer timer;
//      timer.reset();
//      MultiWriter writer("tissue_tumor");

 	  timeloop.execute(onePhasetwoComp);
//      timeloop.executeMultiWriter(onePhasetwoComp, writer);
//      std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

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
  //std::cout << "Please install the UG library." << std::endl;
  std::cout << "Dummy implementation, this test would not compile at the moment." << std::endl;

  return 1;
}
catch (...)
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif
