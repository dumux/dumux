#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
//#include "dumux/material/matrixproperties.hh"
//#include <soilproperties.hh>
#include <soilproperties.hh>
//#include <soilproperties_interface.hh>
#include <dumux/material/twophaserelations.hh>
#include "testproblem_upss.hh"
#include "dumux/timedisc/timeloop.hh"
//#include "dumux/diffusion/fv/fvdiffusion.hh"
#include "dumux/diffusion/fv/fvdiffusionsatupscaled.hh"
//#include "dumux/diffusion/fv/fvdiffusionvelocity.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocityupscaled.hh"
#include "dumux/transport/fv/fvtransport.hh"
//#include "dumux/transport/fv/fvupscaledtransport.hh"
#include "dumux/fractionalflow/impes/impes.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include "dumux/timedisc/expliciteulerstep.hh"
#include "dumux/fractionalflow/variableclass.hh"
#include "../../../dune-subgrid/subgrid/subgrid.hh"
#include "dumux/upscaledsaturation/preprocess/upscalingpreprocess.hh"
//#include "dumux/upscaledsaturation/preprocess/upscalingpreprocess_interface.hh"
//#include "dumux/transport/fv/dispersivecorrection_interface.hh"
#include "dumux/transport/fv/dispersivecorrection.hh"
#include "dumux/transport/fv/convectivecorrection.hh"
//#include "writedispersionmatrix_interface.hh"
#include "writedispersionmatrix.hh"

int main(int argc, char** argv)
{
    try
    {
        // define the problem dimensions
        const int dim = 2;

        std::string arg1(argv[1]);
        std::istringstream is1(arg1);
        double tEnd;
        is1 >> tEnd;

        // create a grid object

        typedef double NumberType;
//        typedef Dune::ALUSimplexGrid<dim,dim> GridType;
        typedef Dune::SGrid<dim,dim> GridType;

//        Dune::GridPtr<GridType> gridpointer("Dxxdata.dgf");
//
//        GridType& grid = *gridpointer;

        typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
        Dune::FieldVector<int,dim> N(10);
        N[0] = 10;
        FieldVector L(0);
        FieldVector H(1);
        H[0] = 1;
        GridType grid(N, L, H);

        grid.globalRefine(2);

        Dune::Water wetmat;
        Dune::Oil nonwetmat;
        //    Dune::HomogeneousSoil<GridType, NumberType> soil;
        //    Dune::HeterogeneousSoil<GridType, NumberType> soil;
        Dune::HeterogeneousSoil<GridType, NumberType> soil(grid, "permeab.dat", false);
        Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wetmat,
                nonwetmat);

        typedef Dune::VariableClass<GridType, NumberType> VariableType;

//        Dune::UpscalingPreprocess<GridType, NumberType> preprocess(grid, wetmat,
//                nonwetmat, soil, materialLaw, 5e6);
        Dune::UpscalingPreprocess<GridType, NumberType> preprocess(grid, wetmat, nonwetmat, soil, materialLaw, 5e4);

        //
//        preprocess.preprocessexecute();

        writeDispersionMatrix(soil);

        double initpress = 0;
        double initsat = 0;
        Dune::FieldVector<double, dim> initvel(0);

        VariableType variables(grid, initsat, initpress, initvel, 0);

        typedef Dune::UpsSProblem<GridType, NumberType, VariableType> Problem;
        Problem problem(variables, wetmat, nonwetmat, soil, materialLaw, L, H,
                false);
//            soil.permeability.vtkout("permeability", grid);
        //
        //    typedef Dune::FVDiffusionVelocity<GridType, NumberType, VariableType> DiffusionType;
        typedef Dune::FVDiffusionVelocityUpscaled<GridType, NumberType, VariableType>
                DiffusionType;
        DiffusionType diffusion(grid, problem, grid.maxLevel());

        typedef Dune::DispersiveCorrection<GridType, NumberType, VariableType>
                Dispersion;
        Dispersion dispersion(problem);

        typedef Dune::ConvectiveCorrection<GridType, NumberType, VariableType>
                ConvectiveCorrection;
        ConvectiveCorrection convectiveCorrection(problem);

        typedef Dune::FVTransport<GridType, NumberType, VariableType>
                TransportType;
//        TransportType transport(grid, problem,0);
        TransportType transport(grid, problem,0,false, dispersion, convectiveCorrection);
        //
        int iterFlag = 2;
        int nIter = 30;
        double maxDefect = 1e-5;
        typedef Dune::IMPES<GridType, DiffusionType, TransportType, VariableType>
                IMPESType;
        IMPESType impes(diffusion, transport, iterFlag, nIter, maxDefect);
        //
        double tStart = 0;
        //    double tEnd = 1.295e8;
        char* fileName("test_upscaledsaturation-testwithcorr-large");
//        char* fileName("test_upscaledsaturation-testnocorr-large");
        int modulo = 1;
        double cFLFactor = 0.95;
        Dune::TimeLoop<GridType, IMPESType> timeloop(tStart, tEnd, fileName,
                modulo, cFLFactor);
        //
        Dune::Timer timer;
        timer.reset();
        timeloop.execute(impes);
        //    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

        return 0;
    } catch (Dune::Exception &e)
    {
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    } catch (...)
    {
        std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
    }
}
