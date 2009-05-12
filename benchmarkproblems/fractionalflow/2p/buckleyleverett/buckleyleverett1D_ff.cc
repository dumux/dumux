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
            coord[i] = coord[i - 1] + (coord[i - 1] - coord[i - 2]) *
            (1 + strfactor);
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

        // IMPES parameters
        int iterFlag = 2;
        int nIter = 30;
        double maxDefect = 1e-5;

        // timeloop parameters
        double tStart = 0;
        //double tEnd = 2.5e9;
        double cFLFactor = 1.0;
        double maxDt = 1e6;
        double firstDt = 1e6;

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

        typedef Dune::VariableClass<GridType, NumberType> VariableType;
        VariableType variables(grid);

        typedef Dune::BuckleyLeverettProblem<GridType, NumberType, VariableType>
                ProblemType;
        ProblemType problem(variables, wPhase, nPhase, soil, materialLaw, Left,
                Right);

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
