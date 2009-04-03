#include "config.h"
#include <iostream>
//#ifdef HAVE_UG
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/timedisc/timeloop.hh"
#include "dumux/fractionalflow/impes/impes_phasepressure.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/transport/fv/fvtransport_wettingphase.hh"
#include "dumux/diffusion/fv/fvdiffusion_pw.hh"
#include "test_decoupled_phasepressure_problem.hh"
#include "test_decoupled_phasepressure_soil.hh"
#include "dumux/fractionalflow/variableclass2p_new.hh"

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions (geometry of problem)
        const int dim=2;
        typedef double NumberType;

        Dune::FieldVector<NumberType, dim> LowerLeft(0);
        Dune::FieldVector<NumberType, dim> UpperRight(300);
        UpperRight[1] = 60;
        Dune::FieldVector<int, dim> elementNum(4);
        elementNum[0] = 30;

        // create a grid object
        typedef Dune::SGrid<dim,dim> GridType;
        GridType grid(elementNum, LowerLeft, UpperRight);

        Dune::gridinfo(grid);

        // time loop parameters
        const double tStart = 0;
        const double tEnd = 4.32e7;
        const double cFLFactor = 0.99;

        // IMPES parameters
        int iterFlag = 2;
        int nIter = 30;
        double maxDefect = 1e-5;
        double omega=0.8;

        // plotting parameters
        const char* fileName = "test_decoupled_phasepressure";
        int modulo = 1;

        // choose fluids and properties
        Dune::Water wPhase;
        Dune::Oil nwPhase;

        typedef Dune::VariableClass<GridType, NumberType> VC;

        VC variables(grid);

        Dune::DecoupledPPTestSoil<GridType, NumberType> soil;

        Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wPhase, nwPhase);

        typedef Dune::DecoupledPPTestProblem<GridType, NumberType, VC> ProblemType;
        ProblemType problem(variables, wPhase, nwPhase, soil, materialLaw, LowerLeft, UpperRight);

        typedef Dune::FVDiffusion<GridType, NumberType, VC, ProblemType> DeprecatedDiffusion;
        DeprecatedDiffusion diffusion(grid, problem);

        typedef Dune::FVTransport<GridType, NumberType, VC, ProblemType> DeprecatedTransport;
        DeprecatedTransport transport(grid, problem);

        typedef Dune::IMPES<GridType, DeprecatedDiffusion, DeprecatedTransport, VC> IMPESType;
        IMPESType impes(diffusion, transport, iterFlag, nIter, maxDefect, omega);

        Dune::TimeLoop<GridType, IMPESType > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);

        Dune::Timer timer;
        timer.reset();
        timeloop.execute(impes);
        std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
