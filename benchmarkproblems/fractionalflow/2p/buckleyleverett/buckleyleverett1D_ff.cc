#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw_deprecated.hh"
#include "dumux/material/brookscoreylaw_deprecated.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"
#include "dumux/transport/fv/fvtransport_deprecated.hh"
#include "dumux/diffusion/fv/fvdiffusion_deprecated.hh"
//#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
//#include "dumux/diffusion/fe/fediffusion.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocity_deprecated.hh"
#include "dumux/fractionalflow/impes/impes_deprecated.hh"
#include "../problemdefinitions/buckleyleveretttransportproblem.hh"
#include "../problemdefinitions/buckleyleverettdiffproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/rungekuttastep.hh"
#include "dumux/fractionalflow/variableclass.hh"
#include "../problemdefinitions/buckleyleverettanalytical.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=1;
    typedef double NumberType;
    typedef GridType::ctype ctype;

    Dune::FieldVector<NumberType, dim> Left(0);
    Dune::FieldVector<NumberType, dim> Right(300);
    if (argc != 2) {
      std::cout << "usage: tEnd" << std::endl;
      return 0;
    }
    std::string arg1(argv[1]);
    std::istringstream is1(arg1);
    double tEnd;
    is1 >> tEnd;

    // create a grid object
    typedef Dune::OneDGrid GridType;

    //deffinition of a stretched grid
//    const int numberofelements = 15;
     const int numberofelements = 30;
//    const int numberofelements = 60;
//    const int numberofelements = 120;

    double strfactor = 0;

    //vector with coordinates
    std::vector<ctype> coord;
    coord.resize(numberofelements+1);
    coord[0]=0;
    coord[1]=1;
    //generate coordinates for a stretched grid
    for (int i=2;i<numberofelements+1;i++){
      coord[i]=coord[i-1]+(coord[i-1]-coord[i-2])*(1+strfactor);
    }

    //scale coordinates to geometry
    for (int i=0;i<numberofelements+1;i++){
      coord[i]*=Right[0]/coord[numberofelements];
      std::cout << "coordinates =  " << coord[i] << std::endl;
    }

    const std::vector<ctype>& coordinates(coord);

    // grid
    GridType grid(coordinates);

    Dune::gridinfo(grid);

    int iterFlag = 2;
    int nIter = 30;
    double maxDefect = 1e-5;

    double tStart = 0;
    //double tEnd = 2.5e9;
    const char* fileName = "buckleyleverett1D";
    int modulo = 1;
    double cFLFactor = 1.0;
    // slope limiter parameters
//    bool reconstruct = true;
//    double alphaMax = 0.8;

    Oil oil(0.2);
    Water water(0.2);
//    Dune::BrooksCoreyLaw materialLaw(water, oil,2,0);
    Dune::LinearLaw materialLaw(water, oil);

    double initpress = 2e5;
    double initsat = 0.2;
    Dune::FieldVector<double,dim> initvel(3e-7);

    typedef Dune::VariableClass<GridType, NumberType> VC;

    VC variables(grid,initsat,initpress,initvel);

    Dune::BLWithAnalytical<GridType, NumberType, VC> transportProblem(variables, materialLaw,Left,Right,cFLFactor);
    Dune::BuckleyLeverettDiffProblem<GridType, NumberType, VC> diffusionProblem(variables, materialLaw,Left,Right);

    typedef Dune::FVTransport<GridType, NumberType, VC> Transport;
    Transport transport(grid, transportProblem, grid.maxLevel());

//    typedef Dune::FVDiffusion<GridType, NumberType, VC> Diffusion;
    typedef Dune::FVDiffusionVelocity<GridType, NumberType, VC> Diffusion;
    Diffusion diffusion(grid, diffusionProblem, grid.maxLevel());


    typedef Dune::IMPES<GridType, Diffusion, Transport, VC> IMPES;
    IMPES fractionalflow(diffusion, transport, iterFlag, nIter, maxDefect);

//    Dune::TimeLoop<GridType, Transport > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);
    Dune::TimeLoop<GridType, IMPES > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);

    Dune::Timer timer;
    timer.reset();
    timeloop.execute(fractionalflow);
//    timeloop.execute(transport);
    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;
    //printvector(std::cout, *fractionalflow, "saturation", "row", 200, 1);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
