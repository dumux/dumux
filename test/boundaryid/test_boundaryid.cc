#include <config.h>
#include <iostream>
#ifdef HAVE_ALBERTA
#include <dune/grid/uggrid.hh>
#include <dune/grid/albertagrid.hh>
#include "gridcheck.cc"
#include "checkgeometryinfather.cc"
#include "checkintersectionit.cc"
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include "dumux/twophase/problems/lensproblem.hh"
#include "lensproblemwithid.hh"
#include "dumux/twophase/fv/boxpwsn.hh"
#include "dumux/timedisc/timeloop.hh"
#include "dumux/material/vangenuchtenlaw.hh"

int main (int argc , char **argv) try
{
    // define the problem dimensions  
    const int dim=2;
    typedef double NumberType; 
    Dune::FieldVector<NumberType, dim> outerLowerLeft(0);
    Dune::FieldVector<NumberType, dim> outerUpperRight(6);
    outerUpperRight[1] = 4;
    Dune::FieldVector<NumberType, dim> innerLowerLeft(1);
    innerLowerLeft[1] = 2;
    Dune::FieldVector<NumberType, dim> innerUpperRight(4);
    innerUpperRight[1] = 3;
    if (argc != 3) {
      std::cout << "usage: test_boundaryid tEnd dt" << std::endl;
      return 0;
    }
    	std::string arg1(argv[1]);
	std::istringstream is1(arg1);
	double tEnd;
	is1 >> tEnd;
	std::string arg2(argv[2]);
	std::istringstream is2(arg2);
	double dt;
	is2 >> dt;


    // create a grid object
	typedef Dune::AlbertaGrid<dim,dim> GridType; 
	//typedef Dune::ALUSimplexGrid<dim,dim> GridType; 

    typedef Dune::BoxPwSn<GridType, NumberType> TwoPhase;
    DNAPL dnapl;
    Water water;
    Dune::VanGenuchtenLaw law(water, dnapl);

    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtrWithID( "grids/rectangle.dgf" );

    // grid reference 
    GridType& gridWithID = *gridPtrWithID;

    Dune::LensProblemWithID<GridType, NumberType> problemWithID(law, outerLowerLeft, outerUpperRight, 
    		innerLowerLeft, innerUpperRight);

    TwoPhase twoPhaseWithID(gridWithID, problemWithID);
    
    Dune::TimeLoop<GridType, TwoPhase> timeloopWithID(0, tEnd, dt, "lenswithid", 1);
    
    timeloopWithID.execute(twoPhaseWithID);
    
    // create grid pointer, GridType is defined by gridtype.hh
    Dune::GridPtr<GridType> gridPtr( "grids/unitcube2.dgf" );

    // grid reference 
    GridType& grid = *gridPtr;

    Dune::LensProblem<GridType, NumberType> problem(law, outerLowerLeft, outerUpperRight, 
    		innerLowerLeft, innerUpperRight);

    TwoPhase twoPhase(grid, problem);
    
    Dune::TimeLoop<GridType, TwoPhase> timeloop(0, tEnd, dt, "lens", 1);
    
    timeloop.execute(twoPhase);

    Dune::BlockVector<Dune::FieldVector<double, 2> > diffVec(*(twoPhase.u));
    diffVec -= *(twoPhaseWithID.u);
    
    //printvector(std::cout, *(twoPhase.u), "without ID", "row", 2, 1, 3);
    //printvector(std::cout, *(twoPhaseWithID.u), "with ID", "row", 2, 1, 3);
    
    std::cout << "difference = " << diffVec.two_norm() << std::endl;
    return 0;
} 
catch (Dune::Exception& e) 
{
    std::cerr << e << std::endl;
    return 1;
 
} 
catch (...) 
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}

#else 

int main (int argc , char **argv) try
{
  std::cout << "Please install the Alberta library." << std::endl;

  return 1;
}
catch (...) 
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif 
