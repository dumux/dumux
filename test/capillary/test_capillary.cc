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
#include "dumux/transport/fv/capillarydiffusion.hh"
#include "dumux/diffusion/fv/fvdiffusion_deprecated.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocity_deprecated.hh"
#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
#include "dumux/fractionalflow/impes/impes_deprecated.hh"
#include "dumux/transport/problems/buckleyleverettproblem.hh"
#include "dumux/transport/problems/simpleproblem.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "dumux/diffusion/problems/heterogeneousproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/rungekuttastep.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;

    // plotting parameters
    const char* fileName = "with_pc";
    int modulo = 3;

    // time loop parameters
    const double tStart = 0;
    const double tEnd = 5e9;
    const double cFLFactor = 0.2;

    // slope limiter parameters
    bool reconstruct = true;
    double alphaMax = 0.8;

    // IMPES parameters
    int iterFlag = 0;
    int nIter = 1;
    double maxDefect = 1e-5;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(1); N[0] = 50;
    FieldVector L(0);
    FieldVector H(300); H[0] = 600;
    GridType grid(N,L,H);

    grid.globalRefine(0);

    Uniform mat;
    //Dune::LinearLaw materialLaw(mat, mat);
    Dune::BrooksCoreyLaw materialLaw(mat, mat,2,10000);
    //Dune::VanGenuchtenLaw materialLaw(mat, mat);

    typedef Dune::VariableClass<GridType, NumberType> VC;

    VC variables(grid);

    Dune::BuckleyLeverettProblem<GridType, NumberType, VC> transportProblem(variables, materialLaw);
    Dune::UniformProblem<GridType, NumberType,VC> diffusionProblem(variables, materialLaw);
    //Dune::HeterogeneousProblem<GridType, NumberType> diffusionProblem(grid, "permeab.dat", false, materialLaw);
    //diffusionProblem.permeability.vtkout("permeability", grid);

    typedef Dune::FVTransport<GridType, NumberType,VC> Transport;
    //Dune::DiffusivePart<GridType, NumberType> diffPart;
    Dune::CapillaryDiffusion<GridType, NumberType,VC> diffPart(diffusionProblem);
    Transport transport(grid, transportProblem, grid.maxLevel(), diffPart, reconstruct, alphaMax);

    typedef Dune::FVDiffusionVelocity<GridType, NumberType,VC> Diffusion;
    Diffusion diffusion(grid, diffusionProblem, grid.maxLevel());

    typedef Dune::IMPES<GridType, Diffusion, Transport,VC> IMPES;
    IMPES fractionalflow(diffusion, transport, iterFlag, nIter, maxDefect);

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
