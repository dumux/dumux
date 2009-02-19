#include "config.h"
#include <iostream>
#define DUMMY
#ifdef DUMMY
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <../../../dune-subgrid/subgrid/subgrid.hh>
#include "dumux/operators/p1operatorextended.hh"
#include <dumux/timedisc/timeloop.hh>
#include <dumux/coupled/boxstokesdarcy.hh>
#include "yxproblem.hh"
#include "../stokes/boxstokes.hh"
#include "boxdiffusion.hh"

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
		  double approximate = solution[globalId];

		  error += (exact - approximate)*(exact - approximate);
	  }

	  return sqrt(error);
}
}

template<int dim>
struct NodeLayout
{
    bool contains(Dune::GeometryType gt) {
        return gt.dim() == 0;
    }
};

int main(int argc, char** argv)
{
  try{
    const int dim=2;
    typedef double NumberType;

    // geometry
    //typedef Dune::ALUSimplexGrid<dim,dim> GridType;
    typedef Dune::SGrid<dim,dim> GridType;
    Dune::GridPtr<GridType> gridPtr( argv[1] );
    GridType& grid = *gridPtr;

    // subdivide grid in subgrids
    typedef Dune::SubGrid<dim,GridType> SubGridType;
    SubGridType subGridStokes(grid);
    SubGridType subGridDarcy(grid);
    subGridStokes.createBegin();
    subGridDarcy.createBegin();
    typedef GridType::Codim<0>::LeafIterator Iterator;
    Iterator eendit = grid.leafend<0>();
    for (Iterator it = grid.leafbegin<0>(); it != eendit; ++it) {
      Dune::GeometryType gt = it->geometry().type();
      const Dune::FieldVector<NumberType,dim>& local = Dune::ReferenceElements<NumberType,dim>::general(gt).position(0, 0);
      Dune::FieldVector<NumberType,dim> global = it->geometry().global(local);
      if (global[1] > 1)
    	  subGridStokes.addPartial(it);
      else
    	  subGridDarcy.addPartial(it);
    }
    subGridStokes.createEnd();
    subGridDarcy.createEnd();

    Dune::YXProblem<SubGridType, NumberType> stokesProblem;
    typedef Dune::LeafP1BoxStokes<SubGridType, NumberType, dim> BoxStokes;
    BoxStokes boxStokes(subGridStokes, stokesProblem);

    DarcyParameters<SubGridType,NumberType> darcyParam;
    typedef Dune::LeafP1BoxDiffusion<SubGridType, NumberType> DarcyModel;
    DarcyModel darcyModel(subGridDarcy, darcyParam);

    typedef Dune::BoxStokesDarcy<BoxStokes,DarcyModel> CoupledModel;
    bool assembleGlobalMatrix = true;
    CoupledModel coupledModel(subGridStokes, boxStokes, subGridDarcy, darcyModel, assembleGlobalMatrix);

    coupledModel.initial();
    coupledModel.assemble();
    printmatrix(std::cout, coupledModel.matrix(), "global stiffness matrix", "row", 11, 4);
    printvector(std::cout, coupledModel.rhs(), "global right hand side", "row", 200, 1, 3);
    coupledModel.solve();
    printvector(std::cout, coupledModel.sol(), "global solution", "row", 200, 1, 3);

    coupledModel.vtkout("test_boxstokesdarcy", 0);

	std::cout << "Darcy discrete error = " << discreteError(subGridDarcy, darcyModel.sol(), darcyParam) << std::endl;

	return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
#else

int main (int argc , char **argv) try
{
  std::cout << "This test is not finished yet." << std::endl;

  return 1;
}
catch (...)
{
    std::cerr << "Generic exception!" << std::endl;
    return 2;
}
#endif
