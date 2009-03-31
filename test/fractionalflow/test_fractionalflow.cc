// $Id$

#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "test_fractionalflow_soilproperties.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include <dumux/material/twophaserelations.hh>
#include "test_fractionalflow_testproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocity.hh"
#include "dumux/transport/fv/fvtransport.hh"
#include "dumux/transport/fv/capillarydiffusion.hh"
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
        Dune::FieldVector<int,dim> N(4); N[0] = 30;
        FieldVector L(0);
        FieldVector H(100 ); H[0] = 300;
        GridType grid(N,L,H);

        grid.globalRefine(0);

        Dune::Water wetmat;
        Dune::Oil nonwetmat;

        Dune::FractionalFlowTestSoil<GridType, NumberType> soil;

        Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wetmat, nonwetmat);

        typedef Dune::VariableClass<GridType, NumberType> VariableType;

        VariableType variables(grid);

        typedef Dune::FractionalFlowTestProblem<GridType, NumberType, VariableType> Problem;
        Problem problem(variables, wetmat, nonwetmat, soil, materialLaw,L, H, false);
        //    soil.permeability.vtkout("permeability", grid);

        typedef Dune::FVDiffusionVelocity<GridType, NumberType, VariableType, Problem> DiffusionType;
        DiffusionType diffusion(grid, problem);

        Dune::CapillaryDiffusion<GridType, NumberType, VariableType, Problem> capillaryDiffusion(problem, soil);

        typedef Dune::FVTransport<GridType, NumberType, VariableType, Problem> TransportType;
        TransportType transport(grid, problem, capillaryDiffusion);

        int iterFlag = 2;
        int nIter = 30;
        double maxDefect = 1e-5;
        typedef Dune::IMPES<GridType, DiffusionType, TransportType, VariableType> IMPESType;
        IMPESType impes(diffusion, transport, iterFlag, nIter, maxDefect);

        double tStart = 0;
        double tEnd = 3e6;
        const char* fileName = "test_fractionalflow";
        int modulo = 10;
        double cFLFactor = 0.95;
        Dune::TimeLoop<GridType, IMPESType > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);

        Dune::Timer timer;
        timer.reset();
        timeloop.execute(impes);
        //    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

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
