// $Id$

#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include "dumux/material/phaseproperties/phaseproperties2p.hh"
#include "dumux/transport/ch/chtransport.hh"
#include "dumux/transport/ch/fractionalw.hh"
#include "dumux/transport/problems/simpleproblem.hh"
#include "dumux/transport/problems/simplenonlinearproblem.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/fractionalflow/variableclass2p.hh"

int main(int argc, char** argv)
{
    try{
        // define the problem dimensions
        const int dim=1;
        typedef double NumberType;

        // time loop parameters
        const double tStart = 0;
        const double tEnd = 1.5e9;
        double dt = 1.5e9;
	double firstDt = dt;
	double maxDt = dt;
        int modulo = 1;

        Dune::FieldVector<NumberType, dim> Left(0);
        Dune::FieldVector<NumberType, dim> Right(600);

        // create a grid object
        typedef Dune::OneDGrid GridType;
        typedef GridType::ctype ctype;
        typedef GridType::LevelGridView GridView;

        //deffinition of a stretched grid
        const int numberofelements = 30;
        double strfactor = 0;

        //vector with coordinates
        std::vector<ctype> coord;
        coord.resize(numberofelements+1);
        coord[0]=0;
        coord[1]=1;

        //generate coordinates for a stretched grid
        for (int i=2;i<numberofelements+1;i++){
            coord[i]=coord[i-1]+(coord[i-1]-coord[i-2])*(1+strfactor);
        }

        //scale coordinates to geometry
        for (int i=0;i<numberofelements+1;i++){
            coord[i]*=Right[0]/coord[numberofelements];
            //std::cout << "coordinates =  " << coord[i] << std::endl;
        }

        const std::vector<ctype>& coordinates(coord);

        // grid
        GridType grid(coordinates);

        grid.globalRefine(0);

        GridView gridView(grid.levelView(0));

        Dune::Uniform mat(0.2);
        Dune::HomogeneousLinearSoil<GridType, NumberType> soil;
        //Dune::HomogeneousNonlinearSoil<GridType, NumberType> soil;
        Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, mat, mat);

        typedef Dune::VariableClass<GridView, NumberType> VC;

        double initsat=0;
        Dune::FieldVector<double,dim>vel(0);
        vel[0] = 1.0/6.0*1e-6;

        VC variables(gridView,initsat,vel);

        Dune::SimpleProblem<GridView, NumberType, VC> problem(variables, mat, mat , soil, materialLaw,Left,Right);
        //Dune::SimpleNonlinearProblem<GridView, NumberType, VC> problem(variables, mat, mat , soil, materialLaw,Left,Right);

        Dune::FractionalW<GridView, NumberType, VC> fractionalW(problem);

        typedef Dune::ChTransport<GridView, NumberType, VC> Transport;

        Transport transport(gridView, problem, fractionalW);

	Dune::RungeKuttaStep<GridType, Transport> timeStep(1);
        Dune::TimeLoop<GridType, Transport > timeloop(tStart, tEnd, dt, "chtransport1D", modulo, maxDt, firstDt, timeStep);

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
