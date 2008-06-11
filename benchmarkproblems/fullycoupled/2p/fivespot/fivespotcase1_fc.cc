#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "../problemdefinitions/fivespotproblem.hh"
#include "dumux/twophase/problems/uniformtwophaseproblem.hh"
#include "dumux/twophase/fv/boxpwsn.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/brookscoreylaw.hh"
#include "dumux/material/vangenuchtenlaw.hh"

int main(int argc, char** argv) {
	try {
		// define the problem dimensions  
		const int dim=2;
		enum {BrooksCorey = 0, VanGenuchten = 1};
		typedef double NumberType;
		Dune::FieldVector<NumberType, dim> LowerLeft(0);
		Dune::FieldVector<NumberType, dim> UpperRight(300);
		UpperRight[1] = 300;
		if (argc < 3) {
			std::cout << "usage: test_twophase tEnd dt" << std::endl;
			return 0;
		}

		if (argc == 4) {
			std::string arg3(argv[3]);
			std::istringstream is3(arg3);
		}
		std::string arg1(argv[1]);
		std::istringstream is1(arg1);
		int tEnd;
		is1 >> tEnd;
		std::string arg2(argv[2]);
		std::istringstream is2(arg2);
		int dt;
		is2 >> dt;
		std::string arg3(argv[3]);
		std::istringstream is3(arg3);
		int bcf;
		is3 >> bcf;

		// create a grid object
		typedef Dune::SGrid<dim,dim> GridType;

		// use unitcube from grids 
		std::stringstream dgfFileName;
		dgfFileName << "grids/unitcube" << GridType :: dimension << ".dgf";

		// create grid pointer, GridType is defined by gridtype.hh
		Dune::GridPtr<GridType> gridPtr( dgfFileName.str() );

		// grid reference 
		GridType& grid = *gridPtr;

		Dune::gridinfo(grid);

		Oil oil(0.2);
		Water water(0.2);
		//    Dune::LinearLaw law (water, oil);
		Dune::BrooksCoreyLaw law(water, oil,2.0,0);
		Dune::FiveSpotCase1Problem<GridType, NumberType> problem(law, LowerLeft, UpperRight, bcf, BrooksCorey);

		typedef Dune::BoxPwSn<GridType, NumberType> TwoPhase;
		TwoPhase twoPhase(grid, problem);

		Dune::TimeLoop<GridType, TwoPhase> timeloop(0, tEnd, dt,
				"fivespotcase1", 1);

		Dune::Timer timer;
		timer.reset();
		timeloop.execute(twoPhase);
		std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

		//printvector(std::cout, *twoPhase.u, "u", "row", 2, 1, 3);

		return 0;
	}
	catch (Dune::Exception &e) {
		std::cerr << "Dune reported error: " << e << std::endl;
	}
	catch (...) {
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}
