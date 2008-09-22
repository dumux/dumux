#include "config.h"
#include <iostream>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
//#include <../../../dune-subgrid/subgrid/subgrid.hh>
#include <dune/disc/stokes/dgstokes.hh>

#include "dumux/operators/p1operatorextended.hh"

#include <dumux/timedisc/timeloop.hh>
#include <dumux/coupled/coupledmodel.hh>
#include "boxdiffusion.hh"

#define DGF

namespace Dune
{
template<int dim>
struct P1Layout
{
	bool contains (Dune::GeometryType gt)
	{
		return gt.dim() == 0;
	}
};
template<class Grid, class Solution, class Problem>
double discreteError(const Grid& grid, const Solution& solution, const Problem& problem)
{
	  enum{dim=Grid::dimension};
		typedef typename Grid::LeafGridView GV;
	    typedef typename GV::IndexSet IS;
	  typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,P1Layout> VM;
		typedef typename GV::template Codim<dim>::Iterator VertexIterator;

	  VM vertexMapper(grid, grid.leafIndexSet());
	  double error = 0.0;
	  const GV& gridview(grid.leafView());

	  VertexIterator endIt = gridview.template end<dim>();
	  VertexIterator it = gridview.template begin<dim>();
	  for (; it != endIt; ++it)
	  {
		  // get exact solution at vertex
		  FieldVector<double,dim> globalCoord = (*it).geometry()[0];
		  double exact = problem.exact(globalCoord);

		  // get approximate solution at vertex
		  int globalId = vertexMapper.map(*it);
		  double approximate = (*solution)[globalId];

		  error += (exact - approximate)*(exact - approximate);
	  }

	  return sqrt(error);
}
}

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;
    typedef double NumberType;
    int refinementSteps = 0;

        Dune::FieldVector<double,dim> length(1.0);
        Dune::FieldVector<int,dim> size(1);
        Dune::FieldVector<bool,dim> periodic(false);
        int overlap = 0;
        typedef Dune::YaspGrid<dim,dim> GridType;
        GridType grid(length, size,periodic,overlap);

    //    Dune::SubGrid subGrid(grid);

    DiffusionParameters<GridType,NumberType> problem;

    typedef Dune::LeafP1BoxDiffusion<GridType, NumberType> Diffusion;
    Diffusion diffusion(grid, problem);

    const int vOrder = 2;
    const int pOrder = 1;
    DGStokesParameters parameters;
    Example<dim, NumberType> exactSolution;
    DirichletBoundary<GridType> dirichletBoundary(exactSolution);
    RightHandSide<GridType> rightHandSide(exactSolution);
    typedef Dune::DGStokes<GridType, vOrder, pOrder> DGStokes;
	DGStokes dGStokes(grid, exactSolution, parameters, dirichletBoundary, rightHandSide, refinementSteps);

    typedef Dune::CoupledModel<Diffusion,DGStokes> CoupledModel;
    CoupledModel coupledModel(grid, diffusion, grid, dGStokes, true);

    coupledModel.initial();
    coupledModel.assemble();
    printmatrix(std::cout, *(diffusion.A), "local stiffness matrix", "row", 11, 4);
    printmatrix(std::cout, coupledModel.matrix(), "global stiffness matrix", "row", 11, 4);

    Dune::TimeLoop<GridType, Diffusion> timeloop(0, 1, 1, "box3d", 1);

    timeloop.execute(diffusion);

    printmatrix(std::cout, *(diffusion.A), "local stiffness matrix", "row", 11, 4);
	std::cout << "discrete error = " << discreteError(grid, *diffusion, problem) << std::endl;
	return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
