#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw.hh"
#include "dumux/material/brookscoreylaw.hh"
#include "dumux/material/vangenuchtenlaw.hh"
#include "dumux/transport/fv/fvtransport.hh"
#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/fractionalflow/impes/impesms.hh"
#include "dumux/transport/problems/buckleyleverettproblem.hh"
#include "dumux/transport/problems/simpleproblem.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "testproblem.hh"
#include "dumux/diffusion/problems/heterogeneousproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/rungekuttastep.hh"
 
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
    
    Dune::BuckleyLeverettProblem<GridType, NumberType> transportProblem(grid, materialLaw);
    Dune::UniformProblem<GridType, NumberType> diffusionProblem(grid, materialLaw);
    //Dune::HeterogeneousProblem<GridType, NumberType> diffusionProblem(grid, "permeab.dat", false, materialLaw);
    //diffusionProblem.permeability.vtkout("permeability", grid);
    //Dune::TestProblem<GridType, NumberType> diffusionProblem(grid, materialLaw);

    typedef Dune::FVTransport<GridType, NumberType> Transport;
    Transport transport(grid, transportProblem, 0);
        
    typedef Dune::FVDiffusion<GridType, NumberType> Diffusion;
    Diffusion diffusion(grid, diffusionProblem, transport.problem, grid.maxLevel());

    int iterFlag = 0; 
    int nIter = 1; 
    double maxDefect = 1e-5;
    typedef Dune::IMPESMS<GridType, Diffusion, Transport> IMPESMS;
    IMPESMS fractionalflow(diffusion, transport, iterFlag, nIter, maxDefect);
    
    double tStart = 0; 
    double tEnd = 1;
    char* fileName("timeloop");
    int modulo = 1; 
    double cFLFactor = 1.0;
    Dune::TimeLoop<GridType, IMPESMS > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);
    
    Dune::Timer timer;
    timer.reset();
    timeloop.execute(fractionalflow);
    printvector(std::cout, *(fractionalflow.diffusion), "pressure", "row", 200, 1);
    printvector(std::cout, *fractionalflow, "saturation", "row", 200, 1);
    printvector(std::cout, fractionalflow.problem.velocity, "velocity", "row", 4, 1);
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
