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
#include "dumux/diffusion/fv/fvdiffusionvelocity_deprecated.hh"
#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
//#include "dumux/fractionalflow/impes/impesms.hh"
#include "dumux/fractionalflow/impes/impes_deprecated.hh"
#include "dumux/transport/problems/buckleyleverettproblem.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "dumux/diffusion/problems/heterogeneousproblem.hh"
#include "dumux/diffusion/problems/levelhetproblem.hh"
//#include "dumux/diffusion/problems/multiscaleproblem.hh"
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
    Dune::FieldVector<int,dim> N(5); N[0] = 5;
    FieldVector L(0);
    FieldVector H(300); H[0] = 300;
    GridType grid(N,L,H);

    grid.globalRefine(2);

    int finelevel = 2;
    int coarselevel = 2;

    Uniform mat;
    Dune::DeprecatedBrooksCoreyLaw materialLaw(mat, mat);

    typedef Dune::VariableClass<GridType, NumberType> VC;

    double initpress = 0;
    double initsat = 0;
    Dune::FieldVector<double, dim> initvel(0);

    VC variables(grid,initpress, initsat,initvel,finelevel);

    Dune::BuckleyLeverettProblem<GridType, NumberType,VC> transportProblem(variables, materialLaw);
    Dune::UniformProblem<GridType, NumberType,VC> diffusionProblem(variables, materialLaw);
//    Dune::LevelHetProblem<GridType, NumberType,VC> diffusionProblem(variables, finelevel, "permeab.dat", false, materialLaw);
//    diffusionProblem.permeability.vtkout("permeability", grid);
//    Dune::MultiscaleProblem<GridType, NumberType , VC> diffProb(grid,diffusionProblem,finelevel,coarselevel, materialLaw);

    typedef Dune::DeprecatedFVTransport<GridType, NumberType,VC> DeprecatedTransport;
    DeprecatedTransport transport(grid, transportProblem, coarselevel);


    typedef Dune::DeprecatedFVDiffusionVelocity<GridType, NumberType, VC> DeprecatedDiffusion;
    DeprecatedDiffusion diffusion(grid, diffusionProblem, finelevel);

//    typedef Dune::MimeticDiffusion<GridType, NumberType> DeprecatedDiffusion;
//    DeprecatedDiffusion diffusion(grid, diffProb, finelevel);

    int iterFlag = 0;
    int nIter = 1;
    double maxDefect = 1e-5;
    typedef Dune::IMPES<GridType, DeprecatedDiffusion, DeprecatedTransport,VC> IMPES;
    IMPES fractionalflow(diffusion, transport, iterFlag, nIter, maxDefect);

//    typedef Dune::IMPESMS<GridType, DeprecatedDiffusion, DeprecatedTransport,VC> IMPESMS;
//    IMPESMS fractionalflow(diffusion, transport, iterFlag, nIter, maxDefect);

    double tStart = 0;
    double tEnd = 2.5e9;
    const char* fileName = "timeloop";
    int modulo = 1;
    double cFLFactor = 0.7;
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
