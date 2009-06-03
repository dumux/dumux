#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "buckleyleverett_soilproperties.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/material/twophaserelations.hh"
#include "dumux/fractionalflow/variableclass2p.hh"
#include "buckleyleverettproblem.hh"
#include "dumux/transport/fv/fvsaturationwetting2p.hh"
#include "dumux/diffusion/fv/fvtotalvelocity2p.hh"
#include "impes_buckleyleverett_analytic.hh"
#include "dumux/timedisc/expliciteulerstep.hh"
#include "dumux/timedisc/timeloop.hh"

int main(int argc, char** argv)
{
    try
    {
        // define the problem dimensions
        const int dim = 1;
        typedef double NumberType;
        typedef GridType::ctype ctype;
        Dune::FieldVector<NumberType, dim> Left(0);
        Dune::FieldVector<NumberType, dim> Right(300);
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
        typedef Dune::OneDGrid GridType;
        typedef GridType::LevelGridView GridView;
        //definition of a stretched grid
        //    const int numberofelements = 15;
        const int numberofelements = 30;
        //    const int numberofelements = 60;
        //    const int numberofelements = 120;
        double strfactor = 0;
        //vector with coordinates
        std::vector<ctype> coord;
        coord.resize(numberofelements + 1);
        coord[0] = 0;
        coord[1] = 1;
        //generate coordinates for a stretched grid
        for (int i = 2; i < numberofelements + 1; i++)
        {
            coord[i] = coord[i - 1] + (coord[i - 1] - coord[i - 2]) * (1
                    + strfactor);
        }
        //scale coordinates to geometry
        for (int i = 0; i < numberofelements + 1; i++)
        {
            coord[i] *= Right[0] / coord[numberofelements];
            std::cout << "coordinates =  " << coord[i] << std::endl;
        }
        const std::vector<ctype>& coordinates(coord);
        // grid
        GridType grid(coordinates);
        Dune::gridinfo(grid);
        GridView gridView(grid.levelView(0));

        // timeloop parameters
        double tStart = 0;
        //double tEnd = 2.5e9;
        double cFLFactor = 1.0;
        double maxDt = 1e6;
        double firstDt = 1e6;

        // IMPES parameters
        int iterFlag = 2;
        int nIter = 30;
        double maxDefect = 1e-5;

        // plotting parameters
        const char* fileName = "buckleyleverett1D";
        int modulo = 1;

        Dune::Water wPhase;
        Dune::Oil nPhase;

        typedef Dune::BuckleyLeverettSoil<GridType, NumberType> SoilType;
        SoilType soil;

        typedef Dune::TwoPhaseRelations<GridType, NumberType>
                TwoPhaseRelationsType;
        TwoPhaseRelationsType materialLaw(soil, wPhase, nPhase);

        typedef Dune::VariableClass<GridView, NumberType> VariableType;
        VariableType variables(gridView);

        typedef Dune::BuckleyLeverettProblem<GridView, NumberType, VariableType>
                ProblemType;
        ProblemType problem(variables, wPhase, nPhase, soil, materialLaw, Left,
                Right);

        typedef Dune::FVTotalVelocity2P<GridView, NumberType, VariableType,
                ProblemType> DiffusionType;
        DiffusionType diffusion(gridView, problem, "pglobal", "Sw");

        typedef Dune::FVSaturationWetting2P<GridView, NumberType, VariableType,
                ProblemType> TransportType;
        TransportType transport(gridView, problem, "vt");

        typedef Dune::IMPESBLAnalytic<GridView, DiffusionType, TransportType,
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
