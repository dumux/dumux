#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
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
#include "dumux/diffusion/mimetic/mimeticdiffusion.hh"
#include "dumux/fractionalflow/impes/impes.hh"
#include "../problemdefinitions/buckleyleveretttransportproblem.hh"
#include "../problemdefinitions/buckleyleverettdiffproblem.hh"
#include "dumux/diffusion/problems/homogeneousproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/rungekuttastep.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv) 
{
   try{
      // define the problem dimensions  
      const int dim=2;
      typedef double NumberType; 
      Dune::FieldVector<NumberType, dim> LowerLeft(0);
      Dune::FieldVector<NumberType, dim> UpperRight(300);
      UpperRight[1]=70;
      if (argc != 2) {
         std::cout << "usage: tEnd" << std::endl;
         return 0;
      }
    	std::string arg1(argv[1]);
	   std::istringstream is1(arg1);
	   int tEnd;
	   is1 >> tEnd;
	  // std::string arg2(argv[2]);
	  // std::istringstream is2(arg2);
	  // int dt;
	  // is2 >> dt;

      // create a grid object
      typedef Dune::SGrid<dim,dim> GridType;

//      // use unitcube from grids 
      std::stringstream dgfFileName;
      dgfFileName << "grids/unitcube" << GridType :: dimension << ".dgf";
//
//      // create grid pointer, GridType is defined by gridtype.hh
      Dune::GridPtr<GridType> gridPtr( dgfFileName.str() );
//
//      // grid reference 
      GridType& grid = *gridPtr;

      Dune::gridinfo(grid);

      Oil oil(0.2);
      Water water(0.2);
      Dune::BrooksCoreyLaw law(water, oil,2,0);
//      Dune::LinearLaw law(water, oil);
      
      typedef Dune::VariableClass<GridType, NumberType> VC;
      
      VC variables(grid);
      
      Dune::BuckleyLeverettTransportProblem<GridType, NumberType, VC> transportProblem(variables, law,LowerLeft,UpperRight);
      Dune::BuckleyLeverettDiffProblem<GridType, NumberType, VC> diffusionProblem(variables, law,LowerLeft,UpperRight);

      typedef Dune::FVTransport<GridType, NumberType, VC> Transport;
      Transport transport(grid, transportProblem, grid.maxLevel());
        
      typedef Dune::FVDiffusionVelocity<GridType, NumberType, VC> Diffusion;
      Diffusion diffusion(grid, diffusionProblem, grid.maxLevel());
//      typedef Dune::MimeticDiffusion<GridType, NumberType> Diffusion;
//          Diffusion diffusion(grid, diffusionProblem, transportProblem);
          
      int iterFlag = 2; 
      int nIter = 30; 
      double maxDefect = 1e-5;
      typedef Dune::IMPES<GridType, Diffusion, Transport, VC> IMPES;
      IMPES fractionalflow(diffusion, transport, iterFlag, nIter, maxDefect);
    
      double tStart = 0; 
      //double tEnd = 2.5e9;
      char* fileName("buckleyleverett");
      int modulo = 1; 
      double cFLFactor = 0.2;
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
