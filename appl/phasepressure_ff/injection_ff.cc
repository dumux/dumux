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
#include "injection_problem.hh"
#include "injection_soil.hh"
#include "dumux/fractionalflow/variableclass2p_new.hh"

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions (geometry of problem)
        const int dim=2;
        typedef double NumberType;
        typedef GridType::ctype ctype;

        Dune::FieldVector<NumberType, dim> lowerLeft(0);
        Dune::FieldVector<NumberType, dim> upperRight(60);
        upperRight[1] = 50;
        Dune::FieldVector<int, dim> elementNum(15);
        elementNum[0] = 20;

        // count number of arguments
        if (argc != 2) {
            std::cout << "usage: test_twophase tEnd" << std::endl;
            return 0;
        }

        // define tEnd
        std::string arg1(argv[1]);
        std::istringstream is1(arg1);
        double tEnd; is1 >> tEnd;

        // create a grid object
        typedef Dune::SGrid<dim,dim> GridType;

        GridType grid(elementNum, lowerLeft, upperRight);

        Dune::gridinfo(grid);

        // time loop parameters
        const double tStart = 0;
        // const double tEnd = 2.5e9;
        const double cFLFactor = 0.99;
        // slope limiter parameters

        // IMPES parameters
        int iterFlag = 2;
        int nIter = 30;
        double maxDefect = 1e-5;
        double omega=0.8;

        // plotting parameters
        const char* fileName = "injection";
        int modulo = 10;

        // choose fluids and properties
        Dune::Water wPhase;
        Dune::Oil nwPhase(860,5e-3);

        typedef Dune::VariableClass<GridType, NumberType> VC;

        VC variables(grid);

        Dune::InjectionSoil<GridType, NumberType> soil;

        Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wPhase, nwPhase);

        typedef Dune::InjectionProblem<GridType, NumberType, VC> ProblemType;
        ProblemType problem(variables, wPhase, nwPhase, soil, materialLaw, lowerLeft, upperRight);

        typedef Dune::FVDiffusion<GridType, NumberType, VC, ProblemType> Diffusion;
        Diffusion diffusion(grid, problem);

        typedef Dune::FVTransport<GridType, NumberType, VC, ProblemType> Transport;
        Transport transport(grid, problem);

        typedef Dune::IMPES<GridType, Diffusion, Transport, VC> IMPES;
        IMPES impes(diffusion, transport, iterFlag, nIter, maxDefect,omega);

        Dune::TimeLoop<GridType, IMPES > timeloop(tStart, tEnd, fileName, modulo, cFLFactor);

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
