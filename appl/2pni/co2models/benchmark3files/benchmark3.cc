//$Id: test_2pni.cc 882 2008-12-04 09:05:55Z melanie $
#include "config.h"
#include <iostream>
//#ifdef HAVE_UG
#include <iomanip>
//#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
//#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
//#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "problem_benchmark3_2pni.hh"
#include "soil_benchmark3_2pni.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/material/matrixproperties.hh"
#include "dumux/material/twophaserelations.hh"
#include "dumux/2pni/fv/boxpwsnte.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/io/vtkmultiwriter.hh"
#include <dune/common/exceptions.hh>

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=3;
        typedef double NumberType;

        // read tEnd and initial time step from console
        if (argc != 4) {
            std::cout << "usage: 2p2cni basefilename tEnd dt" << std::endl;
            return 0;
        }
        std::string arg1(argv[2]);
        std::istringstream is1(arg1);
        double tEnd;
        is1 >> tEnd;
        std::string arg2(argv[3]);
        std::istringstream is2(arg2);
        double dt;
        is2 >> dt;

        // create a grid object
        // typedef Dune::SGrid<dim,dim> GridType;
        //typedef Dune::YaspGrid<dim,dim> GridType;
        //  typedef Dune::UGGrid<dim> GridType;
        //    typedef Dune::ALUSimplexGrid<dim,dim> GridType;
        typedef Dune::ALUCubeGrid<dim,dim> GridType;

        Dune::GridPtr<GridType> gridPointer(argv[1]);
        GridType& grid = *gridPointer;
        //readStarFormat(grid, argv[1]);
        //grid.createLGMGrid(argv[1]);

        Dune::gridinfo(grid);

        // choose fluids
        Dune::Brine wPhase;
        Dune::CO2 nPhase;
        // create soil object
        Dune::SoilBenchmark3_2pni<GridType, NumberType> soil(grid, true, "properties_johansen.dat");
        // create material law object
        Dune::TwoPhaseRelations<GridType, NumberType> law(soil, wPhase, nPhase);

        // create Prolem object
        Dune::ProblemBenchmark3_2pni<GridType, NumberType> problem(wPhase, nPhase, soil, law);

        typedef Dune::VtkMultiWriter<GridType::LeafGridView> MultiWriter;
        typedef Dune::BoxPwSnTe<GridType, NumberType, MultiWriter> TwoPhase;
        TwoPhase twoPhase(grid, problem);

        Dune::TimeLoop<GridType, TwoPhase, true> timeloop(0, tEnd, dt, "dummy", 1);

        Dune::Timer timer;
        timer.reset();
        MultiWriter writer("out-2pni");

        //  for timeloop.executeMultiWriter(twoPhase, writer, true) initial
        //  values are read from restart file data.dgf
        //  at the moment this only works for SGrid in 2D and for ALUCubeGrid in 3D
        timeloop.executeMultiWriter(twoPhase, writer, true);
        std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

        //printvector(std::cout, *twoPhase.u, "u", "row", 2, 1, 3);

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
//#else
//
//int main (int argc , char **argv) try
//{
//  std::cout << "Please install the UG library." << std::endl;
//
//  return 1;
//}

//#endif
