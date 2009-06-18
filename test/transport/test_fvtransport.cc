#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>

#include "dumux/material/fluids/uniform.hh"
#include "dumux/transport/fv/fvsaturationwetting2p.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/fractionalflow/variableclass2p.hh"

#include "simplenonlinearproblem.hh"


int main(int argc, char** argv)
{
	try{
		// define the problem dimensions
		const int dim=2;

		// time loop parameters
		const double tStart = 0;
		const double tEnd = 1.5e9;
		const double cFLFactor = 0.5;
		double maxDT = 1e100;
		int modulo = 10;

		// create a grid object
		typedef double Scalar;
		typedef Dune::SGrid<dim,dim> Grid;
		typedef Grid::LeafGridView GridView;
		Dune::FieldVector<int,dim> N(1); N[0] = 64;
		Dune::FieldVector<Scalar,dim> L(0);
		Dune::FieldVector<Scalar,dim> H(300); H[0] = 600;

		Grid grid(N,L,H);
		GridView gridView(grid.leafView());

		Dune::Uniform fluid;
		Dune::HomogeneousNonlinearSoil<Grid, Scalar> soil;
		Dune::TwoPhaseRelations<Grid, Scalar> materialLaw(soil, fluid, fluid);

		double initsat=0;
		Dune::FieldVector<double,dim> velocity(0);
		velocity[0] = 1.0/6.0*1e-6;
		typedef Dune::VariableClass<GridView, Scalar> VariableClass;
		VariableClass variables(gridView, initsat, velocity);

		Dune::SimpleNonlinearProblem<GridView, Scalar, VariableClass> problem(variables, materialLaw, L, H);

		typedef Dune::FVSaturationWetting2P<GridView, Scalar, VariableClass> Transport;
		Transport transport(gridView, problem, "vt");

		Dune::TimeLoop<Grid, Transport > timeloop(tStart, tEnd, "timeloop", modulo, cFLFactor, maxDT, maxDT);

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
