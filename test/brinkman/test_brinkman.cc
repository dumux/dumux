#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include "dumux/brinkman/fv/fvbrinkman.hh"
#include "brinkmantestproblem.hh"
 
int main(int argc, char** argv) 
{
	try{
		// define the problem dimensions  
		const int dim=3;

		// create a grid object
		typedef double NumberType; 
		typedef Dune::SGrid<dim,dim> GridType; 

		Dune::FieldVector<GridType::ctype,dim> L(0);
		Dune::FieldVector<GridType::ctype,dim> R(1);
		R[0] = 2;
		Dune::FieldVector<int,dim> N(1);
		N[0] = 2;
		GridType grid(N,L,R);

		Dune::BrinkmanTestProblem<GridType, NumberType> problem;

		Dune::FVBrinkman<GridType, NumberType> brinkman(grid, problem);

		printmatrix(std::cout, brinkman.AV, "velocity matrix", "row", 11, 3);
		printvector(std::cout, *brinkman, "pressure", "row", 200, 1, 3);
		brinkman.vtkout("brinkman", 0);

		return 0;
	}
	catch (Dune::Exception &e){
		std::cerr << "Dune reported error: " << e << std::endl;
	}
	catch (...){
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}
