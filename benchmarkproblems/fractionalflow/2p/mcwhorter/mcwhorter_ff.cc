#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw.hh"
#include "dumux/material/brookscoreylaw.hh"
#include "dumux/material/vangenuchtenlaw.hh"
#include "dumux/transport/fv/fvtransport.hh"
#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocity.hh"
//#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
#include "dumux/fractionalflow/impes/impes.hh"
#include "dumux/transport/problems/mcwhortertransportproblem.hh"
#include "dumux/diffusion/problems/mcwhorterdiffproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/rungekuttastep.hh"
#include "dumux/transport/fv/capillarydiffusion.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv) 
{
   try{
      // define the problem dimensions  
      const int dim=2;
      typedef double NumberType; 
     // Dune::FieldVector<NumberType, dim> LowerLeft(0);
     // Dune::FieldVector<NumberType, dim> UpperRight(300);
     // UpperRight[1]=70;
      if (argc != 2) {
         std::cout << "usage: tEnd" << std::endl;
         return 0;
      }
    	std::string arg1(argv[1]);
	   std::istringstream is1(arg1);
	   double tEnd;
	   is1 >> tEnd;
	  
      // create a grid object
      typedef Dune::SGrid<dim,dim> GridType;

      // use unitcube from grids 
      std::stringstream dgfFileName;
      dgfFileName << "grids/unitcube" << GridType :: dimension << ".dgf";

      // create grid pointer, GridType is defined by gridtype.hh
      Dune::GridPtr<GridType> gridPtr( dgfFileName.str() );

      // grid reference 
      GridType& grid = *gridPtr;

      Dune::gridinfo(grid);

      // time loop parameters
      const double tStart = 0;
     // const double tEnd = 2.5e9;
    //  const double cFLFactor = 0.01;
      const double cFLFactor = 0.2;
    
      // slope limiter parameters
      bool reconstruct = true;
      double alphaMax = 0.8;
    
      // IMPES parameters
      int iterFlag = 2; 
      int nIter = 30; 
      double maxDefect = 1e-5;

      // plotting parameters 
      char* fileName("mcwhorter");
      int modulo = 1; 

      Oil oil(0.0);
      Water water(0.0);

      Dune::BrooksCoreyLaw materialLaw(water, oil,2.0,100);
      
      typedef Dune::VariableClass<GridType, NumberType> VC;
      
      VC variables(grid);
    
      Dune::McWhorterTransportProblem<GridType, NumberType, VC> transportProblem(variables, materialLaw);
      Dune::McWhorterDiffProblem<GridType, NumberType, VC> diffusionProblem(variables, materialLaw);

      typedef Dune::FVTransport<GridType, NumberType, VC> Transport;
      Transport transport(grid, transportProblem, grid.maxLevel());
//      Dune::CapillaryDiffusion<GridType, NumberType> diffPart(diffusionProblem);
//      Transport transport(grid, transportProblem, grid.maxLevel(),diffPart,reconstruct, alphaMax, cFLFactor);
        
      typedef Dune::FVDiffusionVelocity<GridType, NumberType, VC> Diffusion;
      Diffusion diffusion(grid, diffusionProblem, grid.maxLevel());

      typedef Dune::IMPES<GridType, Diffusion, Transport, VC> IMPES;
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
