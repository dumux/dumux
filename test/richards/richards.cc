#include "config.h"
#include <iostream>
//#ifdef HAVE_UG
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfalberta.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
//#include "richardsproblemSw.hh"
#include "richardsproblemPw.hh"
//#include "dumux/richards/fv/boxsw.hh"
#include "dumux/richards/fv/boxpw.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/material/twophaserelations.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
	const int dim = 3;
	typedef double NumberType;
    if (argc != 4) {
      std::cout << "usage: richards dgffilename/basefilename tEnd dt" << std::endl;
      return 1;
    }
	std::string arg2(argv[2]);
	std::istringstream is2(arg2);
	double tEnd;
	is2 >> tEnd;
	std::string arg3(argv[3]);
	std::istringstream is3(arg3);
	double dt;
	is3 >> dt;

	//      typedef Dune::ALUSimplexGrid<dim,dim> GridType;
	typedef Dune::ALUCubeGrid<dim,dim> GridType;
	//    typedef Dune::SGrid<dim,dim> GridType;
	    //typedef Dune::ALUSimplexGrid<dim,dim> GridType;
	    //typedef Dune::AlbertaGrid<dim,dim> GridType;
	    //typedef Dune::YaspGrid<dim,dim> GridType;
	    //typedef Dune::UGGrid<dim> GridType;

	// create grid pointer, GridType is defined by gridtype.hh
	Dune::GridPtr<GridType> gridPtr( argv[1] );


    // grid reference
    GridType& grid = *gridPtr;

    Dune::gridinfo(grid);

    Dune::Water wPhase;
    Dune::DNAPL nPhase;
    Dune::RichardsSoil<GridType, NumberType> soil;
    Dune::TwoPhaseRelations<GridType, NumberType> law(soil, wPhase, nPhase);

    Dune::RichardsPwProblem<GridType, NumberType> problem(wPhase, soil, law);

    typedef Dune::BoxPw<GridType, NumberType> Richards;
    Richards richards(grid, problem);

    Dune::TimeLoop<GridType, Richards> timeloop(0, tEnd, dt, "richards", 1);

    Dune::Timer timer;
    timer.reset();
    timeloop.execute(richards);
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
//#else
//
//int main (int argc , char **argv) try
//{
//  std::cout << "Please install the UG library." << std::endl;
//
//  return 1;
//}

//#endif
