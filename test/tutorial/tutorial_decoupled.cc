#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include <tutorial_soilproperties.hh>
#include <dumux/material/twophaserelations.hh>
#include "tutorialproblem_decoupled.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocity.hh"
#include "dumux/transport/fv/fvtransport.hh"
#include "dumux/fractionalflow/impes/impes.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include "dumux/timedisc/expliciteulerstep.hh"
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
    Dune::FieldVector<int,dim> N(10); N[0] = 30;
    FieldVector L(0);
    FieldVector H(300); H[0] = 600;
    GridType grid(N,L,H);

    // define fluid and solid properties and constitutive relationships
    Dune::Water wettingfluid;
    Dune::Oil nonwettingfluid;
    Dune::TutorialSoil<GridType, NumberType> soil;
    Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wettingfluid, nonwettingfluid);

    // create object containing the variables
    typedef Dune::VariableClass<GridType, NumberType> VariableType;
    VariableType variables(grid);

    // create object including the problem definition
    typedef Dune::TutorialProblemDecoupled<GridType, NumberType, VariableType> Problem;
    Problem problem(variables, wettingfluid, nonwettingfluid, soil, materialLaw,L, H);

    // object including the discretisation of the pressure equation
    typedef Dune::FVDiffusionVelocity<GridType, NumberType, VariableType> DiffusionType;
    DiffusionType diffusion(grid, problem, grid.maxLevel());

    // object including the space discretisation of the saturation equation
    typedef Dune::FVTransport<GridType, NumberType, VariableType> TransportType;
    TransportType transport(grid, problem, grid.maxLevel());

    // some parameters used in the IMPES-object
    int iterFlag = 2;
    int nIter = 30;
    double maxDefect = 1e-5;

    // object including the IMPES (IMplicit Pressure Explicit Saturation) Algorithm
    typedef Dune::IMPES<GridType, DiffusionType, TransportType, VariableType> IMPESType;
    IMPESType impes(diffusion, transport, iterFlag, nIter, maxDefect);

    // some parameters needed for the TimeLoop-object
    double tStart = 0;//start simulation at t = tStart
    double tEnd = 2.9e8;//end simulation at t = tEnd
    char* fileName("test_upscaledsaturation");//name of the output files
    int modulo = 1;//define time step interval in which output files are generated
    double cFLFactor = 1;//security factor for the Courant-Friedrichs-Lewy-Criterion

    // create TimeLoop-object
    Dune::TimeLoop<GridType, IMPESType > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);

    Dune::Timer timer;
    timer.reset();

    // start simulation
    timeloop.execute(impes);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
