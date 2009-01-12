// $Id$
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> /*@\label{tutorial-coupled:include-begin}@*/
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "tutorial_soilproperties_coupled.hh"
#include "dumux/material/twophaserelations.hh"
#include "tutorialproblem_coupled.hh"
#include "dumux/twophase/fv/boxpwsn.hh"
#include "dumux/timedisc/timeloop.hh" /*@\label{tutorial-coupled:include-end}@*/
#include "dumux/io/vtkmultiwriter.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2; /*@\label{tutorial-coupled:dim}@*/

    // create a grid object
    typedef double NumberType; /*@\label{tutorial-coupled:grid-begin}@*/
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(10); N[0] = 30;
    FieldVector L(0);
    FieldVector H(300); H[0] = 600;
    GridType grid(N,L,H); /*@\label{tutorial-coupled:grid-end}@*/

    // define fluid and solid properties and constitutive relationships
    Dune::Water wettingfluid; /*@\label{tutorial-coupled:water}@*/
    Dune::Oil nonwettingfluid; /*@\label{tutorial-coupled:oil}@*/
    Dune::TutorialSoil<GridType, NumberType> soil; /*@\label{tutorial-coupled:soil}@*/
    Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wettingfluid, nonwettingfluid);/*@\label{tutorial-coupled:twophaserelations}@*/

    Dune::gridinfo(grid);

    // create object including the problem definition
    typedef Dune::TutorialProblemCoupled<GridType, NumberType> Problem;
    Problem problem(wettingfluid, nonwettingfluid, soil, materialLaw, L, H); /*@\label{tutorial-coupled:problem}@*/

    // create object including the discretisation of the coupled system of
    // equations (oil and water mass balances) with the box method
    typedef Dune::VtkMultiWriter<GridType::LeafGridView> MultiWriter;
    typedef Dune::BoxPwSn<GridType, NumberType, MultiWriter> TwoPhase;
    TwoPhase boxmethod(grid, problem); /*@\label{tutorial-coupled:boxmethod}@*/

    // some parameters needed for the TimeLoop-object
    double tStart = 0; // start simulation at t = tStart
    double tEnd = 1e8; // stop simulation at t = tEnd
//    const char* fileName = "tutorial_coupled"; // name of the output files
//    int modulo = 1; // define time step interval in which output files are generated

    // create TimeLoop-object
    Dune::TimeLoop<GridType, TwoPhase, true> timeloop(tStart, tEnd, 100, "dummy", 1);
//    Dune::TimeLoop<GridType, TwoPhase> timeloop(tStart, tEnd, 100, fileName, modulo); /*@\label{tutorial-coupled:timeloop}@*/

    Dune::Timer timer;
    timer.reset();
    MultiWriter writer("out-tutorial-coupled");
    timeloop.executeMultiWriter(boxmethod, writer); /*@\label{tutorial-coupled:execute}@*/

    // start simulation
//    timeloop.execute(boxmethod); /*@\label{tutorial-coupled:execute}@*/

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
