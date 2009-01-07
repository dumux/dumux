#include "config.h"
#include <iostream>
//#ifdef HAVE_UG
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfalberta.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "co2problem.hh"
#include "dumux/twophase/fv/boxpwsn.hh"
#include "dumux/io/vtkmultiwriter.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/material/twophaserelations.hh"
#include "co2_soilproperties.hh"
#include "dumux/io/vtkmultiwriter.hh"

int main(int argc, char** argv)
{
  try{
	    // define the problem dimensions
	    const int dim=2;
	    typedef double NumberType;
	    double depthBOR = 800.0;  // bottom of reservoir
	    Dune::FieldVector<NumberType, dim> outerLowerLeft(0);
	    Dune::FieldVector<NumberType, dim> outerUpperRight(60);
	    outerUpperRight[1] = 40;

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
	    // typedef Dune::UGGrid<dim> GridType;
		//typedef Dune::ALUSimplexGrid<dim,dim> GridType;

	    Dune::GridPtr<GridType> gridPointer(argv[1]);
	    GridType& grid = *gridPointer;
	    //readStarFormat(grid, argv[1]);
	    //grid.createLGMGrid(argv[1]);

	     Dune::gridinfo(grid);

    Dune::Brine wPhase(1045., 2.535e-4);
    Dune::CO2 nPhase(479., 3.95e-5);
    Dune::Soil<GridType, NumberType> soil;
    Dune::TwoPhaseRelations<GridType, NumberType> law(soil, wPhase, nPhase);

    Dune::CO2Problem<GridType, NumberType> problem(wPhase, nPhase, soil, law, depthBOR);

    typedef Dune::VtkMultiWriter<GridType::LeafGridView> MultiWriter;
    typedef Dune::BoxPwSn<GridType, NumberType, MultiWriter> TwoPhase;
    TwoPhase twoPhase(grid, problem);

//    Dune::TimeLoop<GridType, TwoPhase> timeloop(0, tEnd, dt, "co2-out", 1);
    Dune::TimeLoop<GridType, TwoPhase, true> timeloop(0, tEnd, dt, "dummy", 1);

    Dune::Timer timer;
    timer.reset();
//    timeloop.execute(twoPhase);
    MultiWriter writer("co2-out");

//  for timeloop.executeMultiWriter(twoPhase, writer, true) initial 
//  values are read from restart file data.dgf
//  at the moment this only works for SGrid in 2D and for ALUCubeGrid in 3D    
    timeloop.executeMultiWriter(twoPhase, writer);
    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

    //printvector(std::cout, *twoPhase.u, "u", "row", 2, 1, 3);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
//#else
//
//int main (int argc , char **argv) try
//{
//  std::cout << "Please install the UG library." << std::endl;
//
//  return 1;
//}

//#endif
