//$Id: test_2pni.cc 882 2008-12-04 09:05:55Z melanie $
#include "config.h"
#include <iostream>
//#ifdef HAVE_UG
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfalberta.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "energyproblem.hh"
#include "energysoil.hh"
#include "dumux/2pni/fv/boxpwsnte.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/material/twophaserelations.hh"
#include "dumux/io/vtkmultiwriter.hh"
#include <dune/common/exceptions.hh>

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim = 2;
        typedef double Scalar;

        // read path for grid, total simulation time tEnd and timestep dt
        if (argc != 4)
        {
            std::cout << "usage: 2pni basefilename tEnd dt" << std::endl;
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
        typedef Dune::SGrid<dim,dim> GridType;
        //typedef Dune::YaspGrid<dim,dim> GridType;
        //typedef Dune::UGGrid<dim> GridType;
        //typedef Dune::ALUSimplexGrid<dim,dim> GridType;
        //typedef Dune::ALUCubeGrid<dim,dim> GridType;

        Dune::GridPtr<GridType> gridPointer(argv[1]);
        GridType& grid = *gridPointer;
        //readStarFormat(grid, argv[1]);
        //grid.createLGMGrid(argv[1]);

        Dune::gridinfo(grid);

        // choose fluids
        Dune::Brine wPhase;
        Dune::CO2 nPhase;
        // create soil object
        Dune::TwoPHeatSoil<GridType, Scalar> soil;
        // create material law object
        Dune::TwoPhaseRelations<GridType, Scalar> law(soil, wPhase, nPhase);

        // create Prolem object
        Dune::TwoPHeatProblem<GridType, Scalar> problem(wPhase, nPhase, soil, law);

        typedef Dune::VtkMultiWriter<GridType::LeafGridView> MultiWriter;
        typedef Dune::BoxPwSnTe<GridType, Scalar, MultiWriter> TwoPhase;
        TwoPhase twoPhase(grid, problem);

        Dune::TimeLoop<GridType, TwoPhase, true> timeloop(0, tEnd, dt, "dummy", 1);

        Dune::Timer timer;
        timer.reset();
        MultiWriter writer("out-2pni");

        //  for timeloop.executeMultiWriter(twoPhase, writer, true) restart files are written
        //  for timeloop.executeMultiWriter(twoPhase, writer, true, true, #) initial values are read
        //  from restart file with the number #
        //  at the moment this only works for SGrid in 2D and for ALUCubeGrid in 3D
        timeloop.executeMultiWriter(twoPhase, writer);
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
//#else
//
//int main (int argc , char **argv) try
//{
//  std::cout << "Please install the UG library." << std::endl;
//
//  return 1;
//}

//#endif
