#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw_deprecated.hh"
#include "dumux/material/brookscoreylaw_deprecated.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"
#include "dumux/transport/fv/fvtransport_deprecated.hh"
#include "dumux/transport/fv/capillarydiffusion.hh"
#include "dumux/diffusion/fv/fvdiffusion_deprecated.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocity_deprecated.hh"
//#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
#include "dumux/fractionalflow/impes/impes_deprecated.hh"
#include "../problemdefinitions/mcwhortertransportproblem.hh"
#include "../problemdefinitions/mcwhorteranalytical.hh"
#include "../problemdefinitions/mcwhorterdiffproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/rungekuttastep.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=1;
    typedef double NumberType;
    typedef GridType::ctype ctype;

    Dune::FieldVector<NumberType, dim> Left(0);
    Dune::FieldVector<NumberType, dim> Right(2.6);
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
    const int numberofelements = 56;
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

    // time loop parameters
    const double tStart = 0;
    // const double tEnd = 2.5e9;
    const double cFLFactor = 0.2;
    // slope limiter parameters
    bool reconstruct = true;
    double alphaMax = 0.8;

    // IMPES parameters
    int iterFlag = 2;
    int nIter = 100;
    double maxDefect = 1e-5;
    double omega=1;

    // plotting parameters
    const char* fileName = "mcwhorter1D";
    int modulo = 100;

    Oil oil(0.0);
    Water water(0.0);

    Dune::DeprecatedBrooksCoreyLaw materialLaw(water, oil,2.0,5000);
    //Dune::DeprecatedVanGenuchtenLaw materialLaw(water,oil,3.1257,1.74e-4);
    //Dune::DeprecatedLinearLaw materialLaw(water,oil);

    typedef Dune::VariableClass<GridType, NumberType> VC;

    VC variables(grid);

    Dune::McWWithAnalytical<GridType, NumberType, VC> transportProblem(variables, materialLaw,Left,Right,cFLFactor);
//    Dune::McWhorterTransportProblem<GridType, NumberType, VC> transportProblem(variables, materialLaw,Left,Right);
    Dune::McWhorterDiffProblem<GridType, NumberType, VC> diffusionProblem(variables, materialLaw,Left,Right);

    typedef Dune::DeprecatedFVTransport<GridType, NumberType, VC> DeprecatedTransport;
    //DeprecatedTransport transport(grid, transportProblem, grid.maxLevel());
    Dune::CapillaryDiffusion<GridType, NumberType, VC> diffPart(diffusionProblem);
    DeprecatedTransport transport(grid, transportProblem, grid.maxLevel(),diffPart,reconstruct, alphaMax);

    typedef Dune::DeprecatedFVDiffusionVelocity<GridType, NumberType, VC> DeprecatedDiffusion;
    DeprecatedDiffusion diffusion(grid, diffusionProblem, grid.maxLevel());

    typedef Dune::IMPES<GridType, DeprecatedDiffusion, DeprecatedTransport, VC> IMPES;
    IMPES fractionalflow(diffusion, transport, iterFlag, nIter, maxDefect,omega);

    Dune::TimeLoop<GridType, IMPES > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);

    Dune::Timer timer;
    timer.reset();
    timeloop.execute(fractionalflow);
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
