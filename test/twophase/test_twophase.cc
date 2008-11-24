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
#include "lensproblem.hh"
#include "lenssoil.hh"
#include "dumux/twophase/fv/boxpwsn.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/material/twophaserelations.hh"
#include "dumux/io/vtkmultiwriter.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;
    typedef double NumberType;
    Dune::FieldVector<NumberType, dim> outerLowerLeft(0);
    Dune::FieldVector<NumberType, dim> outerUpperRight(6);
    outerUpperRight[1] = 4;
    Dune::FieldVector<NumberType, dim> innerLowerLeft(1);
    innerLowerLeft[1] = 2;
    Dune::FieldVector<NumberType, dim> innerUpperRight(4);
    innerUpperRight[1] = 3;

    // read tEnd and initial time step from console
    if (argc != 3) {
      std::cout << "usage: test_twophase tEnd dt" << std::endl;
      return 0;
    }
    std::string arg1(argv[1]);
	std::istringstream is1(arg1);
	double tEnd;
	is1 >> tEnd;
	std::string arg2(argv[2]);
	std::istringstream is2(arg2);
	double dt;
	is2 >> dt;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(24); N[1] = 16; // number of cells
    FieldVector L(0); //
    FieldVector H(6); H[1] = 4;
    GridType grid(N,L,H);

    // print some information about the grid
    Dune::gridinfo(grid);

    // choose fluids
    Dune::Water wPhase;
    Dune::DNAPL nPhase;
    // create soil object
    Dune::LensSoil<GridType, NumberType> soil(outerLowerLeft,
    		outerUpperRight, innerLowerLeft, innerUpperRight);
    // create material law object
    Dune::TwoPhaseRelations<GridType, NumberType> law(soil, wPhase, nPhase);

    // create Prolem object
    Dune::LensProblem<GridType, NumberType> problem(wPhase, nPhase, soil, outerLowerLeft,
    		outerUpperRight, innerLowerLeft, innerUpperRight, law);

    typedef Dune::VtkMultiWriter<GridType::LeafGridView> MultiWriter;
    typedef Dune::BoxPwSn<GridType, NumberType, MultiWriter> TwoPhase;
    TwoPhase twoPhase(grid, problem);

    Dune::TimeLoop<GridType, TwoPhase, true> timeloop(0, tEnd, dt, "dummy", 1);

    Dune::Timer timer;
    timer.reset();
    MultiWriter writer("out-twophase");
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
//#else
//
//int main (int argc , char **argv) try
//{
//  std::cout << "Please install the UG library." << std::endl;
//
//  return 1;
//}

//#endif
