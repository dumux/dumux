#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh> // load sgrid definition
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/parabolic/problems/uniformparabolicproblem.hh"
#include "dumux/nonlinear/newtonmethod.hh"
#include "dumux/parabolic/fv/boxnonlinearparabolic.hh"
#include "dumux/parabolic/fe/fenonlinearparabolic.hh"
#include "dumux/timedisc/timeloop.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(1); N[0] = 32;
    FieldVector L(0);
    FieldVector H(300); H[0] = 600;
    GridType grid(N,L,H);

    Dune::UniformParabolicProblem<GridType, NumberType> problem;

    //typedef Dune::FENonlinearParabolic<GridType, NumberType> NonlinearParabolic;
    typedef Dune::BoxNonlinearParabolic<GridType, NumberType> NonlinearParabolic;
    NonlinearParabolic nonlinParabolic(grid, problem);

    Dune::TimeLoop<GridType, NonlinearParabolic> timeloop(0, 5e5, 5e3, "timeloop", 1);

    Dune::Timer timer;
    timer.reset();
    timeloop.execute(nonlinParabolic);
    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

    //printvector(std::cout, *nonlinParabolic.u, "u", "row", 200, 1, 3);

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
