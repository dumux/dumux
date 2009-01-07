#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/common/timer.hh>
#include "dumux/material/linearlaw_deprecated.hh"
#include "dumux/transport/fv/fvtransport_deprecated.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocity_deprecated.hh"
#include "dumux/fractionalflow/impes/impes_deprecated.hh"
#include "dumux/transport/problems/initialballproblem.hh"
#include "dumux/diffusion/problems/gravityproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/fractionalflow/variableclass.hh"
#include "dumux/transport/fv/gravitypart.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;

    // plotting parameters
    const char* fileName = "with_gravity";
    int modulo = 1;

    // time loop parameters
    const double tStart = 0;
    const double tEnd = 3e3;
    const double cFLFactor = 0.99;

    // slope limiter parameters
    bool reconstruct = false;
    double alphaMax = 0.8;

    // IMPES parameters
    int iterFlag = 0;
    int nIter = 1;
    double maxDefect = 1e-5;

    // material law parameters
    //double lambda = 2.0;
    //double p0 = 5e3;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(100); N[0] = 16;
    FieldVector L(0);
    FieldVector H(10); H[0] = 1;
    GridType grid(N,L,H);

    grid.globalRefine(0);

    DNAPL dnapl(0, 2e3, 5e-4);
    Water water(0, 1e3, 5e-4);
    Dune::LinearLaw materialLaw(water, dnapl);

    typedef Dune::VariableClass<GridType, NumberType> VC;
    double initsat = 1;
    VC variables(grid,initsat);

    Dune::InitialBallProblem<GridType, NumberType, VC> transportProblem(variables, grid, materialLaw);
    FieldVector gravity(0);
    gravity[dim-1] = -9.81;
    Dune::GravityProblem<GridType, NumberType, VC> diffusionProblem(variables, &grid, materialLaw, gravity);

    typedef Dune::FVTransport<GridType, NumberType, VC> Transport;
    Dune::GravityPart<GridType, NumberType, VC> gravityPart(diffusionProblem);
    Transport transport(grid, transportProblem, grid.maxLevel(), gravityPart, reconstruct, alphaMax);

    typedef Dune::FVDiffusionVelocity<GridType, NumberType, VC> Diffusion;
    Diffusion diffusion(grid, diffusionProblem, grid.maxLevel());

    typedef Dune::IMPES<GridType, Diffusion, Transport, VC> IMPES;
    IMPES fractionalflow(diffusion, transport, iterFlag, nIter, maxDefect);

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
