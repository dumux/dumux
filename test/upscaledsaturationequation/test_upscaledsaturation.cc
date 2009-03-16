#include "config.h"
#include <iostream>

//#define PARALLEL

#include <dune/common/mpihelper.hh>
//#include <mpi.h>
#include <dune/grid/sgrid.hh> // load sgrid definition
//#include <dune/grid/yaspgrid.hh>
//#include <dune/grid/alugrid.hh>
//#include <dune/grid/io/file/dgfparser/dgfparser.hh>
//#include <dune/grid/io/file/dgfparser/dgfs.hh>
//#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include <soilproperties.hh>
#include <dumux/material/twophaserelations.hh>
#include "testproblem_upss.hh"
#include "dumux/fractionalflow/variableclass.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/diffusion/fv/fvdiffusionsatupscaled.hh"
#include "dumux/diffusion/fv/fvdiffusionvelocityupscaled.hh"
#include "dumux/transport/fv/fvtransport_upscaled.hh"
//include "dumux/transport/fv/fvtransport_upscaled_test.hh"
#include "dumux/fractionalflow/impes/impes_postprocess.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include "../../../dune-subgrid/subgrid/subgrid.hh"
#include "dumux/upscaledsaturation/preprocess/upscalingpreprocess.hh"
#include "dumux/transport/fv/dispersivecorrection.hh"
#include "dumux/transport/fv/convectivecorrection.hh"
#include "writedispersionmatrix.hh"
#include "dumux/timedisc/expliciteulerstep.hh"
#include "dumux/upscaledsaturation/coarsescaleparameters.hh"

int main(int argc, char** argv)
{
    try
    {
#ifdef PARALLEL

        //instantiate MPI
        Dune::MPIHelper::instance(argc, argv);
        Dune::CollectiveCommunication<MPI_Comm> communicator(Dune::MPIHelper::getCommunicator());

        // define the problem dimensions
        const int dim = 2;

        std::string arg1(argv[1]);
        std::istringstream is1(arg1);
        double tEnd;
        is1 >> tEnd;

        // create a grid object

        typedef double Scalar;
        typedef Dune::SGrid<dim,dim> GridType;
        //        typedef Dune::YaspGrid<dim> GridType;

        typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
        Dune::FieldVector<int,dim> N(10);
        N[0] = 10;
        FieldVector L(0);
        FieldVector H(1);
        H[0] = 1;

        GridType grid(N,L,H);

        grid.globalRefine(2);

        Dune::Water wetmat;
        Dune::Oil nonwetmat;

        typedef Dune::HeterogeneousSoil<GridType, Scalar> SoilType;
        SoilType soil(grid, "permeability.dat", false);
        Dune::TwoPhaseRelations<GridType, Scalar> materialLaw(soil, wetmat,
                nonwetmat);
        //                  soil.randomPermeability.vtkout("permeability", grid);

        typedef Dune::VariableClass<GridType, Scalar> VariableType;

        typedef Dune::CoarseScaleParameters<Scalar, Dune::FieldVector<Scalar,dim>, Dune::FieldVector<Scalar,1>, Dune::FieldVector<Scalar,2*dim> >
        CoarseScaleParameters;
        CoarseScaleParameters coarseParameters;

        Dune::UpscalingPreprocess<GridType, Scalar, CoarseScaleParameters>
        preprocess(grid, communicator, wetmat, nonwetmat, soil, coarseParameters,
                materialLaw, 1e6);

        preprocess.preprocessexecute();

        if (communicator.rank() == 0)
        {
            Dune::writeCorrection<CoarseScaleParameters, GridType>(
                    coarseParameters);

            double initpress = 0;
            double initsat = 0;
            Dune::FieldVector<double, dim> initvel(0);

            VariableType variables(grid, initsat, initpress, initvel, 0);

            typedef Dune::UpsSProblem<GridType, Scalar, VariableType,CoarseScaleParameters>
            Problem;
            Problem problem(variables, wetmat, nonwetmat, soil,
                    coarseParameters, materialLaw, L, H, false);
            //          soil.randomPermeability.vtkout("permeability", grid);
            //
            //    typedef Dune::FVDiffusionVelocity<GridType, Scalar, VariableType> DiffusionType;
            typedef Dune::FVDiffusionVelocityUpscaled<GridType, Scalar, VariableType,Problem>
            DiffusionType;
            DiffusionType diffusion(grid, problem);

            typedef Dune::DispersiveCorrection<GridType, Scalar, VariableType, Problem>
            DispersiveCorrection;
            DispersiveCorrection dispersiveCorrection(problem);

            typedef Dune::ConvectiveCorrection<GridType, Scalar, VariableType, Problem>
            ConvectiveCorrection;
            ConvectiveCorrection convectiveCorrection(problem);

            typedef Dune::FVTransport<GridType, Scalar, VariableType,Problem>
            TransportType;
            //                  TransportType transport(grid, problem,0,false);
            TransportType transport(grid, problem, 0, false,
                    dispersiveCorrection, convectiveCorrection);
            //
            int iterFlag = 2;
            int nIter = 30;
            int modulo = 1;
            double maxDefect = 1e-5;
            typedef Dune::IMPESPostProcess<GridType, DiffusionType, TransportType, VariableType>
            IMPESType;
            IMPESType impes(diffusion, transport, modulo, iterFlag, nIter,
                    maxDefect);
            //
            double tStart = 0;
            //    double tEnd = 1.295e8;
            char* fileName("test_upscaledsaturation-testwithcorr-large");
            //      char* fileName("test_upscaledsaturation-testnocorr-large");

            double cFLFactor = 0.8;
            double maxDt = 1e100;
            double firstDt = 1e100;
            Dune::ExplicitEulerStep<GridType, IMPESType> timestep;
            Dune::TimeLoop<GridType, IMPESType> timeloop(tStart, tEnd, fileName,
                    modulo, cFLFactor, maxDt, firstDt, timestep);
            //
            Dune::Timer timer;
            timer.reset();
            timeloop.execute(impes);
        }
#else
        // define the problem dimensions
        const int dim = 2;

        std::string arg1(argv[1]);
        std::istringstream is1(arg1);
        double tEnd;
        is1 >> tEnd;

        // create a grid object

        typedef double Scalar;
        typedef Dune::SGrid<dim,dim> GridType;

        typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
        Dune::FieldVector<int,dim> N(10);
        N[0] = 10;
        FieldVector L(0);
        FieldVector H(1);
        H[0] = 1;
        GridType grid(N, L, H);

        grid.globalRefine(3);

        Dune::Water wetmat;
        Dune::Oil nonwetmat;

        typedef Dune::HeterogeneousSoil<GridType, Scalar> SoilType;
        SoilType soil(grid, "permeability.dat", false);
        Dune::TwoPhaseRelations<GridType, Scalar> materialLaw(soil, wetmat,
                nonwetmat);

        typedef Dune::VariableClass<GridType, Scalar> VariableType;

        typedef Dune::CoarseScaleParameters<Scalar, Dune::FieldVector<Scalar,dim>, Dune::FieldVector<Scalar,1>, Dune::FieldVector<Scalar,2*dim> >
                CoarseScaleParameters;
        CoarseScaleParameters coarseParameters;

        Dune::UpscalingPreprocess<GridType, Scalar, CoarseScaleParameters>
                preprocess(grid, wetmat, nonwetmat, soil,
                        coarseParameters, materialLaw, 1e6);

//        preprocess.preprocessexecute();

//        Dune::writeCorrection<CoarseScaleParameters, GridType>(coarseParameters);

        double initpress = 0;
        double initsat = 0;
        Dune::FieldVector<double, dim> initvel(0);

        VariableType variables(grid, initsat, initpress, initvel, 0);

        typedef Dune::UpsSProblem<GridType, Scalar, VariableType,CoarseScaleParameters>
                Problem;
        Problem problem(variables, wetmat, nonwetmat, soil, coarseParameters,
                materialLaw, L, H, false);
        //          soil.randomPermeability.vtkout("permeability", grid);
        //
        //    typedef Dune::FVDiffusionVelocity<GridType, Scalar, VariableType> DiffusionType;
        typedef Dune::FVDiffusionVelocityUpscaled<GridType, Scalar, VariableType,Problem>
                DiffusionType;
        DiffusionType diffusion(grid, problem);

        typedef Dune::DispersiveCorrection<GridType, Scalar, VariableType, Problem>
                DispersiveCorrection;
        DispersiveCorrection dispersiveCorrection(problem);

        typedef Dune::ConvectiveCorrection<GridType, Scalar, VariableType, Problem>
                ConvectiveCorrection;
        ConvectiveCorrection convectiveCorrection(problem);

        typedef Dune::FVTransport<GridType, Scalar, VariableType,Problem>
                TransportType;
                          TransportType transport(grid, problem,0,false);
//        TransportType transport(grid, problem, 0, false, dispersiveCorrection,
//                convectiveCorrection);
        //
        int iterFlag = 2;
        int nIter = 30;
        int modulo = 1;
        double maxDefect = 1e-5;
        typedef Dune::IMPESPostProcess<GridType, DiffusionType, TransportType, VariableType>
                IMPESType;
        IMPESType impes(diffusion, transport, modulo, iterFlag, nIter,
                maxDefect);
        //
        double tStart = 0;
        //    double tEnd = 1.295e8;
//        char* fileName("test_upscaledsaturation-testwithcorr-large");
              char* fileName("test_upscaledsaturation-testnocorr-large");

        double cFLFactor = 0.8;
        double maxDt = 1e100;
        double firstDt = 1e100;
        Dune::ExplicitEulerStep<GridType, IMPESType> timestep;
        Dune::TimeLoop<GridType, IMPESType> timeloop(tStart, tEnd, fileName,
                modulo, cFLFactor, maxDt, firstDt, timestep);
        //
        Dune::Timer timer;
        timer.reset();
        timeloop.execute(impes);

#endif
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
