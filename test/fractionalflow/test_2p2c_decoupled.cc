#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> 
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw_deprecated.hh"
#include "dumux/material/brookscoreylaw_deprecated.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"
//#include "dumux/transport/fv/fvtransport2p2c.hh"
#include "dumux/diffusion/fv/fvdiffusion.hh"
//#include "dumux/fractionalflow/impes/impesms.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/expliciteulerstep.hh"
//#include "dumux/transport/problems/testproblem_2p2c.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "dumux/diffusion/problems/levelhetproblem.hh"
//#include "dumux/diffusion/problems/fourspotproblem.hh" 
#include "dumux/transport/fv/efendiev_dispersion.hh"

int main(int argc, char** argv) 
{
  try{
//    // define the problem dimensions  
//    const int dim=2;
//
//    // create a grid object
//    typedef double NumberType; 
//    typedef Dune::SGrid<dim,dim> GridType; 
//    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector; 
//    Dune::FieldVector<int,dim> N(5); N[0] = 5;                   
//    FieldVector L(0); 
//    FieldVector H(300); H[0] = 300; 
//    GridType grid(N,L,H);  
// 
//    grid.globalRefine(3);
//
//    Water wetmat;
//    Air nonwetmat;
//    Dune::LinearLaw materialLaw(wetmat, nonwetmat);
//    
//    HenryWaterAir henry;
//    
//
////  Dune::UniformProblem<GridType, NumberType> diffusionProblem(grid, materialLaw);
//    typedef Dune::LevelHetProblem<GridType, NumberType> DiffProb;
////  typedef Dune::FourSpotProblem<GridType, NumberType> DiffProb;
//	DiffProb diffusionProblem(grid, 3, "permeab.dat", true, materialLaw);
//    diffusionProblem.permeability.vtkout("permeability", grid);
//    
//    typedef Dune::Testproblem_2p2c<GridType, NumberType> TransProb;
//    TransProb transportProblem(henry, grid, materialLaw);
//    
//    Dune::DiffusivePart<GridType,NumberType> diffusivePart;
//        
//    typedef Dune::FVDiffusion<GridType, NumberType> Diffusion;
//    Diffusion diffusion(grid, diffusionProblem, transportProblem, grid.maxLevel());
//
//    transportProblem.diffusion = &diffusion;
//    
//    typedef Dune::FVTransport2p2c<GridType, NumberType> Transport;
//    Transport transport(grid, transportProblem, grid.maxLevel(), diffusivePart, false);
//    
//    int iterFlag = 0; 
//    int nIter = 1; 
//    double maxDefect = 1e-5;
//    typedef Dune::IMPESMS<GridType, Diffusion, Transport> IMPESMS;
//    IMPESMS fractionalflow(diffusion, transport, iterFlag, nIter, maxDefect);
//    
//    double tStart = 0; 
//    double tEnd = 1e5;
//    char* fileName("timeloop");
//    int modulo = 1; 
//    double cFLFactor = 0.99;
//    
//    Dune::ExplicitEulerStep<GridType, IMPESMS> timestep;
//    
//    transport.initialguess();
//    //diffusion.pressure(transport.sat);
//    diffusion.pressure();
//    
//    Dune::TimeLoop<GridType, IMPESMS > timeloop(tStart, tEnd, fileName, modulo, cFLFactor, 1e100, 1e100, timestep);
//    
//    Dune::Timer timer;
//    timer.reset();
//    timeloop.execute(fractionalflow);
////    printvector(std::cout, *(fractionalflow.diffusion), "pressure", "row", 200, 1);
////    printvector(std::cout, *fractionalflow, "saturation", "row", 200, 1);
////    printvector(std::cout, fractionalflow.problem.velocity, "velocity", "row", 4, 1);
//    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;
    
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
