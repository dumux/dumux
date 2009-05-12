#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "../problemdefinitions/buckleyleverett_soilproperties.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/material/twophaserelations.hh"
#include "../problemdefinitions/buckleyleverettproblem.hh"
#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocity.hh"
#include "dumux/transport/fv/fvtransport.hh"
#include "dumux/fractionalflow/impes/impes_buckleyleverett_analytic.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/expliciteulerstep.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv)
{
    try
    {
        // define the problem dimensions
        const int dim = 2;
        typedef double NumberType;
        Dune::FieldVector<NumberType, dim> LowerLeft(0);
        Dune::FieldVector<NumberType, dim> UpperRight(300);
        UpperRight[1] = 70;
        if (argc != 2)
        {
            std::cout << "usage: tEnd" << std::endl;
            return 0;
        }
        std::string arg1(argv[1]);
        std::istringstream is1(arg1);
        double tEnd;
        is1 >> tEnd;
        // std::string arg2(argv[2]);
        // std::istringstream is2(arg2);
        // int dt;
        // is2 >> dt;

        // create a grid object
        typedef Dune::SGrid<dim, dim> GridType;
        // use unitcube from grids
        std::stringstream dgfFileName;
        dgfFileName << "grids/unitcube" << GridType::dimension << ".dgf";
        // create grid pointer, GridType is defined by gridtype.hh
        Dune::GridPtr<GridType> gridPtr(dgfFileName.str());
        // grid reference
        GridType& grid = *gridPtr;
        Dune::gridinfo(grid);

        // timeloop parameters
        double tStart = 0;
        //double tEnd = 2.5e9;
        double cFLFactor = 1;
        double maxDt = 1e6;
        double firstDt = 10;

        // IMPES parameters
        int iterFlag = 2;
        int nIter = 30;
        double maxDefect = 1e-5;

        // plotting parameters
        const char* fileName = "buckleyleverett";
        int modulo = 1;

        Dune::Water wPhase;
        Dune::Oil nPhase;

        typedef Dune::BuckleyLeverettSoil<GridType, NumberType> SoilType;
        SoilType soil;

        typedef Dune::TwoPhaseRelations<GridType, NumberType>
                TwoPhaseRelationsType;
        TwoPhaseRelationsType materialLaw(soil, wPhase, nPhase);

        typedef Dune::VariableClass<GridType, NumberType> VariableType;
        VariableType variables(grid);

        typedef Dune::BuckleyLeverettProblem<GridType, NumberType, VariableType>
                ProblemType;
        ProblemType problem(variables, wPhase, nPhase, soil, materialLaw,
                LowerLeft, UpperRight);

        typedef Dune::FVDiffusionVelocity<GridType, NumberType, VariableType,
                ProblemType> DiffusionType;
        DiffusionType diffusion(grid, problem);

        typedef Dune::FVTransport<GridType, NumberType, VariableType,
                ProblemType> TransportType;
        TransportType transport(grid, problem);

        typedef Dune::IMPESBLAnalytic<GridType, DiffusionType, TransportType,
                VariableType> IMPESType;
        IMPESType impes(diffusion, transport, iterFlag, nIter, maxDefect);

        Dune::ExplicitEulerStep<GridType, IMPESType> timestep;
        Dune::TimeLoop<GridType, IMPESType> timeloop(tStart, tEnd, fileName,
                modulo, cFLFactor, maxDt, firstDt, timestep);
        Dune::Timer timer;
        timer.reset();
        timeloop.execute(impes);
        std::cout << "timeloop.execute took " << timer.elapsed() << " seconds"
                << std::endl;
        //printvector(std::cout, *fractionalflow, "saturation", "row", 200, 1);

        return 0;
    } catch (Dune::Exception &e)
    {
        std::cerr << "Dune reported error: " << e << std::endl;
    } catch (...)
    {
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
