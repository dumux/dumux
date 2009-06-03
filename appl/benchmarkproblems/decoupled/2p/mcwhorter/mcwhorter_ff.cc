#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "mcwhorter_soilproperties.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/material/twophaserelations.hh"
#include "dumux/fractionalflow/variableclass2p.hh"
#include "mcwhorterproblem.hh"
#include "dumux/transport/fv/fvsaturationwetting2p.hh"
#include "dumux/transport/fv/capillarydiffusion.hh"
#include "dumux/diffusion/fv/fvwettingvelocity2p.hh"
#include "impes_mcwhorter_analytic.hh"
#include "dumux/timedisc/expliciteulerstep.hh"
#include "dumux/timedisc/timeloop.hh"

int main(int argc, char** argv)
{
    try
    {
        // define the problem dimensions
        const int dim = 2;
        typedef double NumberType;
        Dune::FieldVector<NumberType, dim> LowerLeft(0);
        Dune::FieldVector<NumberType, dim> UpperRight(2.6);
        UpperRight[1] = 1;
        Dune::FieldVector<int, dim> cellNumbers(4);
        cellNumbers[0] = 26;
        if (argc != 2)
        {
            std::cout << "usage: tEnd" << std::endl;
            return 0;
        }
        std::string arg1(argv[1]);
        std::istringstream is1(arg1);
        double tEnd;
        is1 >> tEnd;

        // create a grid object
        typedef Dune::SGrid<dim, dim> GridType;
        typedef GridType::LevelGridView GridView;

        // grid reference
        GridType grid(cellNumbers,LowerLeft,UpperRight);
        Dune::gridinfo(grid);

        GridView gridView(grid.levelView(0));

        // time loop parameters
        const double tStart = 0;
        // const double tEnd = 2.5e9;
        // const double cFLFactor = 0.01;
        const double cFLFactor = 0.1;
        double maxDt = 1e100;
        double firstDt = 1e100;

        // IMPES parameters
        int iterFlag = 2;
        int nIter = 30;
        double maxDefect = 1e-5;

        // plotting parameters
        const char* fileName = "mcwhorter";
        int modulo = 1;

        Dune::Water wPhase(1000,0.001);
        Dune::Oil nPhase(1000,0.001);

        typedef Dune::McWhorterSoil<GridType, NumberType> SoilType;
        SoilType soil;

        Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wPhase,
                nPhase);

        typedef Dune::VariableClass<GridView, NumberType> VariableType;
        VariableType variables(gridView);

        typedef Dune::McWhorterProblem<GridView, NumberType, VariableType>
                ProblemType;
        ProblemType problem(variables, wPhase, nPhase, soil, materialLaw,
                UpperRight);

        typedef Dune::FVSaturationWetting2P<GridView, NumberType, VariableType,
                ProblemType> Transport;
        Dune::CapillaryDiffusion<GridView, NumberType, VariableType,
                ProblemType> diffPart(problem, soil);
        Transport transport(gridView, problem, "vw");

        typedef Dune::FVWettingPhaseVelocity2P<GridView, NumberType, VariableType,
                ProblemType> Diffusion;
        Diffusion diffusion(gridView, problem, "pn", "Sw");

        typedef Dune::IMPESMcWAnalytic<GridView, Diffusion, Transport,
                VariableType> IMPESType;
//        typedef Dune::IMPES<GridView, Diffusion, Transport,
//                        VariableType> IMPESType;
        IMPESType impes(diffusion, transport, iterFlag, nIter, maxDefect);

        Dune::ExplicitEulerStep<GridType, IMPESType> timestep;
        Dune::TimeLoop<GridType, IMPESType> timeloop(tStart, tEnd, fileName,
                modulo, cFLFactor, maxDt, firstDt, timestep);
        Dune::Timer timer;
        timer.reset();
        timeloop.execute(impes);
        std::cout << "timeloop.execute took " << timer.elapsed() << " seconds"
                << std::endl;
        // printvector(std::cout, *fractionalflow, "saturation", "row", 200, 1);

        return 0;
    } catch (Dune::Exception &e)
    {
        std::cerr << "Dune reported error: " << e << std::endl;
    } catch (...)
    {
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
