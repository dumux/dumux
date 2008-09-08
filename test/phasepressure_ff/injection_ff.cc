#include "config.h"
#include <iostream>
//#ifdef HAVE_UG
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"
#include "dumux/material/brookscoreylaw_deprecated.hh"
#include "dumux/fractionalflow/impes/impes_deprecated.hh"
#include "dumux/material/properties.hh"
#include "dumux/transport/fv/fvtransport_phase.hh"
#include "dumux/transport/fv/capillarydiffusion.hh"
#include "dumux/diffusion/fv/fvdiffusion_pw.hh"
#include "dumux/diffusion/fv/fvdiffusionwettingphasevelocity_pw.hh"
#include "injectiontransportproblem.hh"
#include "injectiondiffproblem.hh"
#include "dumux/timedisc/rungekuttastep.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions (geometry of problem)
    const int dim=2;
    typedef double NumberType;
    typedef GridType::ctype ctype;

    Dune::FieldVector<NumberType, dim> outerLowerLeft(0);
    Dune::FieldVector<NumberType, dim> outerUpperRight(60);
    outerUpperRight[1] = 50;
    double depthBOR = 1000.0;

    // for defining e.g. a lense
    Dune::FieldVector<NumberType, dim> innerLowerLeft(4);
    innerLowerLeft[1] = 0.0;
    Dune::FieldVector<NumberType, dim> innerUpperRight(6);
    innerUpperRight[1] = 0.5;

    // count number of arguments
    if (argc != 2) {
      std::cout << "usage: test_twophase tEnd" << std::endl;
      return 0;
    }

    // define tEnd
    std::string arg1(argv[1]);
    std::istringstream is1(arg1);
    double tEnd; is1 >> tEnd;

    // create a grid object
    typedef Dune::SGrid<dim,dim> GridType;
    //typedef Dune::YaspGrid<dim,dim> GridType;
//    typedef Dune::UGGrid<dim> GridType;

    // use unitcube from grids (UGGrid)
    std::stringstream dgfFileName;
    dgfFileName << "grids/unitcube"
//    dgfFileName << "grids/unitcube"
    	<< GridType :: dimension << ".dgf";

    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr( dgfFileName.str() );

    // grid reference
    GridType& grid = *gridPtr;

    Dune::gridinfo(grid);

    // time loop parameters
    const double tStart = 0;
    // const double tEnd = 2.5e9;
    const double cFLFactor = 0.01;
    // slope limiter parameters
    bool reconstruct = false;
    double alphaMax = 0.8;

    // IMPES parameters
    int iterFlag = 2;
    int nIter = 500;
    double maxDefect = 1e-5;
    double omega=0.8;

    // plotting parameters
    char* fileName("injection");
    int modulo = 20;

    // choose fluids and properties
     Water wPhase(0,1000); CO2 nwPhase(0,630,6e-5);

//    Dune::LinearLaw materialLaw(wphase,nwphase,10000);
    Dune::BrooksCoreyLaw materialLaw(wPhase, nwPhase,2,10000);

    typedef Dune::VariableClass<GridType, NumberType> VC;

    VC variables(grid);

    Dune::InjectionTransportProblem<GridType, NumberType, VC> transportProblem(variables, materialLaw, outerLowerLeft, outerUpperRight,
    		innerLowerLeft, innerUpperRight);
    Dune::InjectionDiffProblem<GridType, NumberType, VC> diffusionProblem(variables, materialLaw,outerLowerLeft, outerUpperRight,
    		innerLowerLeft, innerUpperRight,depthBOR, true);

    typedef Dune::FVTransport<GridType, NumberType, VC> Transport;
    Transport transport(grid, transportProblem, grid.maxLevel(),reconstruct, alphaMax);

    typedef Dune::FVDiffusionVelocity<GridType, NumberType, VC> Diffusion;
    Diffusion diffusion(grid, diffusionProblem,  grid.maxLevel());

    typedef Dune::IMPES<GridType, Diffusion, Transport, VC> IMPES;
    IMPES fractionalflow(diffusion, transport, iterFlag, nIter, maxDefect,omega);

    Dune::TimeLoop<GridType, IMPES > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);

    Dune::Timer timer;
    timer.reset();
    timeloop.execute(fractionalflow);
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
//#else
//
//int main (int argc , char **argv) try
//{
//  std::cout << "Please install the UG library." << std::endl;
//
//  return 1;
//}
//catch (...)
//{
//    std::cerr << "Generic exception!" << std::endl;
//    return 2;
//}
//#endif
