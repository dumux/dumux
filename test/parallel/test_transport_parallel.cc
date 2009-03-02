#include "config.h"
#include <iostream>
#include <iomanip>
#undef DUMMY
#ifdef DUMMY
#include <dune/common/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include "dumux/material/properties.hh"
#include "dumux/material/linearlaw_deprecated.hh"
#include "dumux/material/brookscoreylaw_deprecated.hh"
#include "dumux/material/vangenuchtenlaw_deprecated.hh"
#include "dumux/transport/fv/fvtransport_deprecated.hh"
#include "dumux/transport/fv/capillarydiffusion.hh"
#include "dumux/transport/problems/buckleyleverettproblem.hh"
#include "dumux/transport/problems/simpleproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/timedisc/rungekuttastep.hh"
#include "dumux/diffusion/problems/uniformproblem.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv)
{
    try{
        Dune::MPIHelper::instance(argc, argv);

        // define the problem dimensions
        const int dim=2;

        // time loop parameters
        const double tStart = 0;
        const double tEnd = 2.5e9;
        const double cFLFactor = 1.0;
        double maxDT = 1e100;
        int modulo = 10;

        // create a grid object
        typedef double NumberType;
        Dune::FieldVector<double,dim> lowerLeft(0);
        Dune::FieldVector<double,dim> length(300); length[0] = 600;
        Dune::FieldVector<int,dim> size(1); size[0] = 64;
        Dune::FieldVector<bool,dim> periodic(false);
        int overlap = 1;
        typedef Dune::YaspGrid<dim,dim> GridType;
        GridType grid(MPI_COMM_WORLD, length, size, periodic, overlap);

        Uniform mat(0.2);
        //Dune::VanGenuchtenLaw materialLaw(mat, mat);
        //Dune::BrooksCoreyLaw materialLaw(mat, mat);
        Dune::LinearLaw materialLaw(mat, mat);

        typedef Dune::VariableClass<GridType, NumberType> VC;

        double initsat=0;
        double initpress=0;
        Dune::FieldVector<double,dim>vel(0);
        vel[0] = 1.0/6.0*1e-6;

        VC variables(grid,initsat,initpress,vel);

        Dune::SimpleProblem<GridType, NumberType, VC> problem(variables, materialLaw,lowerLeft,length);

        typedef Dune::FVTransport<GridType, NumberType, VC> Transport;

        Transport transport(grid, problem, grid.maxLevel());


        Dune::TimeLoop<GridType, Transport > timeloop(tStart, tEnd, "transport_parallel", modulo, cFLFactor, maxDT, maxDT);

        Dune::Timer timer;
        timer.reset();
        timeloop.execute(transport);
        double elapsedTime = timer.elapsed();
        elapsedTime = grid.comm().max(elapsedTime);
        if (grid.comm().rank() == 0)
            std::cout << "timeloop.execute took " << elapsedTime << " seconds" << std::endl;

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
#else

int main (int argc , char **argv) try
 {
     std::cout << "This test is not finished yet." << std::endl;
     //  std::cout << "Please install MPI." << std::endl;

     return 1;
 }
 catch (...)
 {
     std::cerr << "Generic exception!" << std::endl;
     return 2;
 }
#endif
