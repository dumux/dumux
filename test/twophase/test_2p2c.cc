#include "config.h"
#include <iostream>

#define DUMMY
//#undef DUMMY
#ifdef DUMMY
//#ifdef HAVE_UG

#include <iomanip>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
//#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "injectionproblem.hh"
#include "dumux/2p2c/fv/box2p2c.hh"
#include "dumux/timedisc/timeloop.hh"

#include "dumux/material/phaseproperties/phaseproperties_waterair.hh"
#include "dumux/material/matrixproperties.hh"
#include "dumux/material/twophaserelations.hh"

#include "dumux/material/multicomponentrelations.hh"
#include "dumux/io/vtkmultiwriter.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions (geometry of problem)
    const int dim=2;
    typedef double NumberType;
    Dune::FieldVector<NumberType, dim> outerLowerLeft(0.0);
    Dune::FieldVector<NumberType, dim> outerUpperRight(60.0);
    outerUpperRight[1] = 40.0;
    double depthBOR = 1000.0;

    // for defining e.g. a lense
    Dune::FieldVector<NumberType, dim> innerLowerLeft(4);
    innerLowerLeft[1] = 0.0;
    Dune::FieldVector<NumberType, dim> innerUpperRight(6);
    innerUpperRight[1] = 0.5;

    if (argc != 4) {
      std::cout << "usage: 2p2cni grid tEnd dt" << std::endl;
      return 0;
    }
    // define tEnd
    std::string arg1(argv[2]);
	std::istringstream is1(arg1);
	double tEnd;
	is1 >> tEnd;
    // define dt
	std::string arg2(argv[3]);
	std::istringstream is2(arg2);
	double dt;
	is2 >> dt;

    // create a grid object
    //typedef Dune::SGrid<dim,dim> GridType;
    //typedef Dune::YaspGrid<dim,dim> GridType;
    typedef Dune::UGGrid<dim> GridType;

    // use grid defined in the arguments
    Dune::GridPtr<GridType> gridPointer(argv[1]);
    // grid reference
    GridType& grid = *gridPointer;

    Dune::gridinfo(grid);

    // choose fluids and properties
    Dune::Liq_WaterAir wPhase;
    Dune::Gas_WaterAir nPhase;
    Dune::Injectionsoil<GridType, NumberType> soil;

    Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wPhase, nPhase);
    Dune::CWaterAir multicomp(wPhase, nPhase);

    // create problem properties and geometry
    Dune::InjectionProblem<GridType, NumberType> problem(wPhase, nPhase, soil, outerLowerLeft,
    		outerUpperRight, innerLowerLeft, innerUpperRight, depthBOR, materialLaw, multicomp);

    // create two-phase two-component problem
    typedef Dune::VtkMultiWriter<GridType> MultiWriter;
    typedef Dune::Box2P2C<GridType, NumberType, MultiWriter> TwoPhaseTwoComp;
    TwoPhaseTwoComp twoPhasetwoComp(grid, problem);

    Dune::TimeLoop<GridType, TwoPhaseTwoComp> timeloop(0, tEnd, dt, "lens", 1);

    Dune::Timer timer;
    timer.reset();
    Dune::VtkMultiWriter<GridType> writer("out2p2c");
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
//  std::cout << "Please install the UG library." << std::endl;
  std::cout << "Dummy implementation, this test would not compile at the moment." << std::endl;

  return 1;
}
catch (...)
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif
