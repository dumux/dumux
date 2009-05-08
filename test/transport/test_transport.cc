#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/transport/fv/fvtransport.hh"
#include "dumux/transport/problems/simpleproblem.hh"
#include "dumux/transport/problems/simplenonlinearproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/fractionalflow/variableclass.hh"

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=2;

        // time loop parameters
        const double tStart = 0;
        const double tEnd = 2.5e9;
        const double cFLFactor = 0.5;
        double maxDT = 1e100;
        int modulo = 1;

        // slope limiter parameters
        //    bool reconstruct = true;
        //    double alphaMax = 0.8;

        // create a grid object
        typedef double NumberType;
        typedef Dune::SGrid<dim,dim> GridType;
        typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
        Dune::FieldVector<int,dim> N(1); N[0] = 64;
        FieldVector L(0);
        FieldVector H(300); H[0] = 600;
        GridType grid(N,L,H);

        std::stringstream dgfFileName;
        dgfFileName << "grids/unitcube" << GridType :: dimension << ".dgf";

        grid.globalRefine(0);

        Dune::Uniform mat(0.2);
        //Dune::HomogeneousLinearSoil<GridType, NumberType> soil;
        Dune::HomogeneousNonlinearSoil<GridType, NumberType> soil;
        Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, mat, mat);

        typedef Dune::VariableClass<GridType, NumberType> VC;

        double initsat=0;
        double initpress=0;
        Dune::FieldVector<double,dim>vel(0);
        vel[0] = 1.0/6.0*1e-6;

        VC variables(grid,initsat,initpress,vel);

        //Dune::SimpleProblem<GridType, NumberType, VC> problem(variables, mat, mat , soil, materialLaw,L,H);
        Dune::SimpleNonlinearProblem<GridType, NumberType, VC> problem(variables, mat, mat , soil, materialLaw,L,H);

        typedef Dune::TransportProblem<GridType, NumberType, VC> TransportProblem;
        typedef Dune::FVTransport<GridType, NumberType, VC, TransportProblem> Transport;

        Transport transport(grid, problem);


        Dune::TimeLoop<GridType, Transport > timeloop(tStart, tEnd, "timeloop", modulo, cFLFactor, maxDT, maxDT);

        timeloop.execute(transport);

        printvector(std::cout, variables.saturation, "saturation", "row", 200, 1);

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
