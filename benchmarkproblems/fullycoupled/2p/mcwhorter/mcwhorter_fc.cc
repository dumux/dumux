# include "config.h"
#include <iostream>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "../problemdefinitions/mcwhorterproblem.hh"
#include "dumux/twophase/problems/uniformtwophaseproblem.hh"
#include "dumux/twophase/fv/boxpnsw.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/brookscoreylaw_deprecated.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"
#include "dumux/io/vtkmultiwriter.hh"
#include"../problemdefinitions/mcwhorteranalytical.hh"

int main(int argc, char** argv)
{
   try{
      // define the problem dimensions
      const int dim=2;
      typedef double NumberType;
      Dune::FieldVector<NumberType, dim> LowerLeft(0);
      Dune::FieldVector<NumberType, dim> UpperRight(2.6);
      UpperRight[1] = 1;
      if (argc < 3) {
         std::cout << "usage: test_twophase tEnd dt" << std::endl;
         return 0;
      }

      std::string arg1(argv[1]);
       std::istringstream is1(arg1);
       int tEnd;
       is1 >> tEnd;
       std::string arg2(argv[2]);
       std::istringstream is2(arg2);
       int dt;
       is2 >> dt;

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


      Oil oil(0);
      Water water(0);
      Uniform mat;
      Dune::BrooksCoreyLaw law(water, oil,2,5000);
//      Dune::LinearLaw law(water, oil);
//      Dune::VanGenuchtenLaw law(water, oil);

      //calculate with analytical solution
      Dune::McWWithAnalytical<GridType, NumberType> problem(grid,law, LowerLeft, UpperRight);

      //calculate without analytical solution
//      Dune::McWhorterProblem<GridType, NumberType> problem(law, LowerLeft, UpperRight);

      typedef Dune::VtkMultiWriter<GridType::LeafGridView> MultiWriter;
      typedef Dune::BoxPnSw<GridType, NumberType, MultiWriter> TwoPhase;
      TwoPhase twoPhase(grid, problem);

      Dune::TimeLoop<GridType, TwoPhase, true> timeloop(0, tEnd, dt,
      "mcwhorter", 1);

      Dune::Timer timer;
      timer.reset();
      
      MultiWriter writer("mcworther_fc");
      timeloop.executeMultiWriter(twoPhase, writer);
      std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

      //printvector(std::cout, *twoPhase.u, "u", "row", 2, 1, 3);

      return 0;
    }
   catch (Dune::Exception &e){
       std::cerr << "Dune reported error: " << e << std::endl;
   }
   catch (...){
      std::cerr << "Unknown exception thrown!" << std::endl;
   }
}
