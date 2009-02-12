#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "lensproblem.hh"
#include "lenssoilwithid.hh"
//#include "dumux/twophase/fv/boxpwsn.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/material/twophaserelations.hh"
#include "dumux/io/vtkmultiwriter.hh"
#include "dumux/new_models/2p/pwsnboxmodel.hh"

int main(int argc, char** argv)
{
    try{
        typedef double Scalar;

        if (argc != 4) {
            std::cout << "usage: test_elementid basefilename tEnd dt" << std::endl;
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

        // instantiate and run the concrete problem (defined in
        // lensproblem.hh)
        Dune::Lens::PwSnLensProblem<Scalar> problem(dt, tEnd, argv[1]);
        if (!problem.simulate())
            return 2;

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
