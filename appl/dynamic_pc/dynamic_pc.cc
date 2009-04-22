#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/material/matrixproperties.hh"
#include "dumux/transport/fv/fvtransport_dynamic_pc.hh"
#include "dynamic_pc_problem.hh"
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
        const double cFLFactor = 1.0;
        double maxDT = 1e100;
        int modulo = 1;

        // create a grid object
        typedef double NumberType;
        typedef Dune::SGrid<dim,dim> GridType;
        typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
        Dune::FieldVector<int,dim> N(1); N[0] = 64;
        FieldVector L(0);
        FieldVector H(300); H[0] = 600;
        GridType grid(N,L,H);

        grid.globalRefine(0);

        Dune::Uniform mat;
        Dune::HomogeneousSoil<GridType, NumberType> soil;
        Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, mat, mat);

        typedef Dune::VariableClass<GridType, NumberType> VariableType;

        double initsat=0;
        double initpress=0;
        Dune::FieldVector<double,dim>vel(0);
        vel[0] = 1.0/6.0*1e-6;

        VariableType variables(grid,initsat,initpress,vel);

        typedef Dune::DynamicPcProblem<GridType, NumberType, VariableType> Problem;
        Problem problem(variables, mat, mat , soil, materialLaw,L,H);

        typedef Dune::TransportProblem<GridType, NumberType, VariableType> TransportProblem;
        typedef Dune::FVTransport<GridType, NumberType, VariableType, TransportProblem> TransportType;

//        Dune::CapillaryDiffusion<GridType, NumberType, VariableType, Problem> capillaryDiffusion(problem, soil);

        TransportType transport(grid, problem);


        Dune::TimeLoop<GridType, TransportType > timeloop(tStart, tEnd, "timeloop", modulo, cFLFactor, maxDT, maxDT);

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
