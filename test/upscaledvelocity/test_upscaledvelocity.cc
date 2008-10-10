#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw_deprecated.hh"
#include "dumux/material/brookscoreylaw_deprecated.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"
#include "dumux/transport/fv/fvtransport_deprecated.hh"
#include "dumux/diffusion/fv/fvdiffusion_deprecated.hh"
#include "dumux/fractionalflow/impes/impesms.hh"
#include "dumux/transport/problems/buckleyleverettproblem.hh"
#include "dumux/transport/problems/simpleproblem.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
//#include "testproblem.hh"
#include "dumux/diffusion/problems/heterogeneousproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/rungekuttastep.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(1); N[0] = 2;
    FieldVector L(0);
    FieldVector H(300); H[0] = 600;
    GridType grid(N,L,H);

    grid.globalRefine(1);

    Uniform mat;
    Dune::LinearLaw materialLaw(mat, mat);

    typedef Dune::VariableClass<GridType, NumberType> VC;
    double initsat=0;
    double initpress=0;
    Dune::FieldVector<double,dim>vel(0);
    vel[0] = 1.0/6.0*1e-6;
    VC variables(grid,initsat,initpress,vel);

    Dune::BuckleyLeverettProblem<GridType, NumberType, VC> transportProblem(variables, materialLaw, grid.maxLevel());
    Dune::UniformProblem<GridType, NumberType, VC> diffusionProblem(variables, materialLaw, true);
    //Dune::HeterogeneousProblem<GridType, NumberType> diffusionProblem(grid, "permeab.dat", false, materialLaw);
    //diffusionProblem.permeability.vtkout("permeability", grid);
    //Dune::TestProblem<GridType, NumberType> diffusionProblem(grid, materialLaw);

    typedef Dune::FVTransport<GridType, NumberType, VC> Transport;
    Transport transport(grid, transportProblem, 0);

    typedef Dune::FVDiffusion<GridType, NumberType, VC> Diffusion;
    Diffusion diffusion(grid, diffusionProblem, grid.maxLevel());

    int iterFlag = 0;
    int nIter = 1;
    double maxDefect = 1e-5;
    typedef Dune::IMPESMS<GridType, Diffusion, Transport, VC> IMPESMS;
    IMPESMS fractionalflow(diffusion, transport, iterFlag, nIter, maxDefect);

    double tStart = 0;
    double tEnd = 1;
    const char* fileName = "timeloop";
    int modulo = 1;
    double cFLFactor = 1.0;
    Dune::TimeLoop<GridType, IMPESMS > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);

    Dune::Timer timer;
    timer.reset();
    timeloop.execute(fractionalflow);
    printvector(std::cout, fractionalflow.diffproblem.variables.pressure, "pressure", "row", 200, 1);
    printvector(std::cout, fractionalflow.transproblem.variables.saturation, "saturation", "row", 200, 1);
    printvector(std::cout, fractionalflow.diffproblem.variables.velocity, "velocity", "row", 4, 1);
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
