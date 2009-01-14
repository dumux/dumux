//$Id$

#define DUNE_DEVEL_MODE
#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/phaseproperties/phaseproperties_waterair.hh"
#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>
#include "dumux/transport/problems/testproblem2p2cni.hh"
#include "dumux/transport/fv/decoupled2p2cni.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include "dumux/timedisc/expliciteulerstep.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(20); N[0] = 20;
    FieldVector L(0);
    FieldVector H(300); H[0] = 300;
    GridType grid(N,L,H);

    grid.globalRefine(0);

    double tStart = 0;
    double tEnd = 2.5e5;
    int modulo = 1;
    double cFLFactor = 0.9;

    Dune::Liq_WaterAir wetmat;
    Dune::Gas_WaterAir nonwetmat;
    Dune::HomogeneousSoil<GridType, NumberType> soil;
//    Dune::HeterogeneousSoil<GridType, NumberType> soil(grid, "permeab.dat", true);
    Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wetmat, nonwetmat);

    Dune::VariableClass2p2cni<GridType,NumberType> var(grid);

    typedef Dune::TestProblem2p2cni<GridType, NumberType> TransProb;
    TransProb problem(var, wetmat, nonwetmat, soil, materialLaw, false);

    Dune::DiffusivePart<GridType, NumberType> diffPart;
    const Dune::Upwind<NumberType> numFl;

    typedef Dune::Decoupled2p2cni<GridType, NumberType> ModelType;
    ModelType model(grid, problem, grid.maxLevel(), diffPart, false, 0.8, numFl, "CG");

    Dune::ExplicitEulerStep<GridType, ModelType> timestep;
        Dune::TimeLoop<GridType, ModelType > timeloop(tStart, tEnd, "2p2cni", modulo, cFLFactor, 1e100, 1e100, timestep);

        Dune::Timer timer;
        timer.reset();
        timeloop.execute(model);
        std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

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
