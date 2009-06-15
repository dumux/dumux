// $Id$

#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/yaspgrid.hh> // load yaspgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include "dumux/material/fluids/uniform.hh"
#include "dumux/transport/ch/chtransport.hh"
#include "dumux/transport/ch/fractionalw.hh"
#include "simplenonlinearproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/fractionalflow/variableclass2p.hh"

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=2;

        // time loop parameters
        const double tStart = 0;
        const double tEnd = 1.5e9;
        double dt = 3e8;
	double firstDt = dt;
	double maxDt = dt;
        int modulo = 1;
	double slLengthFactor = 10;

        // create a grid object
        typedef double NumberType;

        Dune::FieldVector<int,dim> N(1); N[0] = 100;
        Dune::FieldVector<NumberType,dim> L(0);
        Dune::FieldVector<NumberType,dim> H(300); H[0] = 600;


        typedef Dune::SGrid<dim,dim> GridType;
        typedef GridType::LevelGridView GridView;
        GridType grid(N,L,H);
        //typedef Dune::YaspGrid<dim> GridType;
        //GridType grid(H, N, Dune::FieldVector<bool,dim>(false), 0);

        grid.globalRefine(0);

        GridView gridView(grid.levelView(0));

        Dune::Uniform mat(0.2);
        //Dune::HomogeneousLinearSoil<GridType, NumberType> soil;
        Dune::HomogeneousNonlinearSoil<GridType, NumberType> soil;
        Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, mat, mat);

        typedef Dune::VariableClass<GridView, NumberType> VC;

        double initsat=0;
        Dune::FieldVector<double,dim>vel(0);
        vel[0] = 1.0/6.0*1e-6;

        VC variables(gridView,initsat,vel);

        //Dune::SimpleProblem<GridType, NumberType, VC> problem(variables, mat, mat , soil, materialLaw,L,H);
        Dune::SimpleNonlinearProblem<GridView, NumberType, VC> problem(variables, mat, mat , soil, materialLaw,L,H);

        Dune::FractionalW<GridView, NumberType, VC> fractionalW(problem);

        typedef Dune::ChTransport<GridView, NumberType, VC> Transport;

        Transport transport(gridView, problem, fractionalW, slLengthFactor);

	Dune::RungeKuttaStep<GridType, Transport> timeStep(1);
        Dune::TimeLoop<GridType, Transport > timeloop(tStart, tEnd, dt, "chtransport", modulo, maxDt, firstDt, timeStep);

        timeloop.execute(transport);

        printvector(std::cout, variables.saturation(), "saturation", "row", 200, 1);

        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
