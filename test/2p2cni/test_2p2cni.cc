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
#include "waterairproblem.hh"

#include "dumux/material/phaseproperties/phaseproperties_waterair.hh"
#include "dumux/material/matrixproperties.hh"
#include "dumux/material/twophaserelations.hh"

#include "dumux/2p2cni/fv/box2p2cni.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/io/readstarformat.cc"
#include "dumux/io/vtkmultiwriter.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;
    typedef double NumberType;
    double depthBOR = 1000.0;  // bottom of reservoir
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
    //typedef Dune::UGGrid<dim> GridType;
	//typedef Dune::ALUSimplexGrid<dim,dim> GridType;

    Dune::GridPtr<GridType> gridPointer(argv[1]);
    GridType& grid = *gridPointer;
    //readStarFormat(grid, argv[1]);
    //grid.createLGMGrid(argv[1]);

     Dune::gridinfo(grid);

     // choose fluids and properties
     Dune::Liq_WaterAir wPhase;
     Dune::Gas_WaterAir nPhase;
     Dune::HomogeneousSoil<GridType, NumberType> soil;

     Dune::TwoPhaseRelations<GridType, NumberType>
     materialLaw(soil, wPhase, nPhase);

     Dune::CWaterAir multicomp(wPhase, nPhase);

     Dune::WaterAirProblem<GridType, NumberType> problem(wPhase, nPhase, soil,
    		 materialLaw, multicomp, depthBOR);

     typedef Dune::VtkMultiWriter<GridType::LeafGridView> MultiWriter;
     typedef Dune::Box2P2CNI<GridType, NumberType, MultiWriter> TwoPhase;
     TwoPhase twoPhase(grid, problem);

     Dune::TimeLoop<GridType, TwoPhase> timeloop(0, tEnd, dt, "out2p2cni", 1);

     Dune::Timer timer;
     timer.reset();
     MultiWriter writer("out2p2cni");

	 //timeloop.execute(twoPhase);
     timeloop.executeMultiWriter(twoPhase, writer);
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
