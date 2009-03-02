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
#include "dmtproblem.hh"
#include "dmtsoil.hh"
#include "dumux/twophase/fv/boxpwsn.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/material/twophaserelations.hh"
#include "dumux/io/vtkmultiwriter.hh"

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=2;
        typedef double NumberType;

        // read tEnd and initial time step from console
        if (argc != 4) {
            std::cout << "usage: test_dmt basefilename tEnd dt" << std::endl;
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
        typedef Dune::UGGrid<dim> GridType;
        //typedef Dune::ALUSimplexGrid<dim,dim> GridType;

        Dune::GridPtr<GridType> gridPtr(argv[1]);

        GridType& grid = *gridPtr;

        // choose fluids
        Dune::Methane wPhase(0, 1.02e-5); // variable density and constant viscosity
        Dune::Water nPhase(999.67, 1.3e-3); // constant density and viscosity
        // create soil object
        Dune::DMTSoil<GridType, NumberType> soil(gridPtr);
        // create material law object
        Dune::TwoPhaseRelations<GridType, NumberType> law(soil, wPhase, nPhase);

        // create Prolem object
        Dune::DMTProblem<GridType, NumberType> problem(wPhase, nPhase, soil, law);

        typedef Dune::VtkMultiWriter<GridType::LeafGridView> MultiWriter;
        typedef Dune::BoxPwSn<GridType, NumberType, MultiWriter> TwoPhase;
        TwoPhase twoPhase(grid, problem);

        Dune::TimeLoop<GridType, TwoPhase, true> timeloop(0, tEnd, dt, "dummy", 1);

        Dune::Timer timer;
        timer.reset();
        MultiWriter writer("out-dmt");
        timeloop.executeMultiWriter(twoPhase, writer);
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
