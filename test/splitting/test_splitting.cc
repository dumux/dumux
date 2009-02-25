#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw_deprecated.hh"
#include "dumux/material/brookscoreylaw_deprecated.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"
#include "dumux/transport/fv/fvtransport_deprecated.hh"
#include "dumux/transport/fv/fvsplittedtransport.hh"
#include "dumux/transport/fv/capillarydiffusion.hh"
#include "dumux/transport/problems/buckleyleverettproblem.hh"
//#include "dumux/transport/problems/simpleproblem.hh"
#include "dumux/transport/problems/simpleparabolicproblem.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/rungekuttastep.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;

    // time loop parameters
    const double tStart = 0;
    const double tEnd = 2.5e9;
    double cFLFactor = 0.3;
    int modulo = 10;

    // slope limiter parameters
    bool reconstruct = false;
    double alphaMax = 0.8;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(1); N[0] = 64;
    FieldVector L(0);
    FieldVector H(300); H[0] = 600;
    GridType grid(N,L,H);

    std::stringstream dgfFileName;
    dgfFileName << "grids/unitcube" << GridType :: dimension << ".dgf";

    grid.globalRefine(0);

    Uniform mat;
    //Dune::DeprecatedVanGenuchtenLaw materialLaw(mat, mat);
    Dune::DeprecatedBrooksCoreyLaw materialLaw(mat, mat);
    //Dune::DeprecatedLinearLaw materialLaw(mat, mat);

    typedef Dune::VariableClass<GridType, NumberType> VC;
    double initsat=0;
    double initpress=0;
    Dune::FieldVector<double,dim>vel(0);
    vel[0] = 1.0/6.0*1e-6;
    VC variables(grid,initsat,initpress,vel);

    typedef Dune::DeprecatedFVTransport<GridType, NumberType, VC> HyperbolicPart;
    typedef Dune::DeprecatedFVTransport<GridType, NumberType, VC> ParabolicPart;
    typedef Dune::FVSplittedTransport<GridType, NumberType, VC> SplittedTransport;

    Dune::DiffusivePart<GridType, NumberType> diffPart;
    //Dune::SimpleProblem<GridType, NumberType, VC> hyperbolicProblem(variables, materialLaw, true);
    Dune::SimpleParabolicProblem<GridType, NumberType, VC> hyperbolicProblem(variables, materialLaw, true);
    HyperbolicPart hyperbolicPart(grid, hyperbolicProblem, grid.maxLevel(), diffPart, reconstruct, alphaMax);

    Dune::UniformProblem<GridType, NumberType, VC> diffusionProblem(variables, materialLaw, true);
    //Dune::DiffusivePart<GridType, NumberType> diffPart2;
    Dune::CapillaryDiffusion<GridType, NumberType, VC> diffPart2(diffusionProblem); // use CapillaryDiffusion
    Dune::SimpleParabolicProblem<GridType, NumberType, VC> parabolicProblem(variables, grid, materialLaw);
    ParabolicPart parabolicPart(grid, parabolicProblem, grid.maxLevel(), diffPart2);

    SplittedTransport splittedTransport(grid, hyperbolicPart, parabolicPart);

    Dune::TimeLoop<GridType, SplittedTransport > timeloop(tStart, tEnd, "splitted", modulo, cFLFactor);

    timeloop.execute(splittedTransport);

    printvector(std::cout, *splittedTransport, "saturation", "row", 200, 1);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
