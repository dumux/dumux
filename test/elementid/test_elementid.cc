#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "lensproblem.hh"
#include "lenssoilwithid.hh"
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

    if (argc != 4) {
      std::cout << "usage: test_elementid basefilename tEnd dt" << std::endl;
      return 0;
    }
        std::string arg1(argv[2]);
    std::istringstream is1(arg1);
    double tEnd;
    is1 >> tEnd;
    std::string arg2(argv[3]);
    std::istringstream is2(arg2);
    double dt;
    is2 >> dt;

    typedef Dune::UGGrid<dim> GridType;
  // create grid pointer, GridType is defined by gridtype.hh
  Dune::GridPtr<GridType> gridPtr( argv[1] );

  // grid reference
  GridType& grid = *gridPtr;
  grid.globalRefine(2);

    // print some information about the grid
    Dune::gridinfo(grid);

    // choose fluids
    Dune::Water wPhase;
    Dune::DNAPL nPhase;
    // create soil object
    Dune::LensSoilWithId<GridType, NumberType> soil(gridPtr, outerLowerLeft,
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
    MultiWriter writer("out-elementid");
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
