#include "config.h"
#include <iostream>
#define DUMMY
#ifdef DUMMY
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <../../../dune-subgrid/subgrid/subgrid.hh>

#include "dumux/operators/p1operatorextended.hh"

#include <dumux/timedisc/timeloop.hh>
#include <dumux/coupled/coupleddiffusion.hh>
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
          FieldVector<double,dim> globalCoord = (*it).geometry().corner(0);
          double exact = problem.exact(globalCoord);

          // get approximate solution at vertex
          int globalId = vertexMapper.map(*it);
          double approximate = solution[globalId];

          error = std::max(fabs(exact - approximate), error);
      }

      return (error);
}
}

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;
    typedef double NumberType;

    Dune::FieldVector<double,dim> length(1.0);
    length[0] = 2.0;
    Dune::FieldVector<int,dim> size(32);
    size[0] = 2*size[1];
    Dune::FieldVector<bool,dim> periodic(false);
    int overlap = 0;
    typedef Dune::YaspGrid<dim,dim> GridType;
    GridType grid(length,size,periodic,overlap);

    typedef Dune::SubGrid<dim,GridType> SubGridType;
    SubGridType subGridLeft(grid);
    SubGridType subGridRight(grid);
    subGridLeft.createBegin();
    subGridRight.createBegin();
    typedef GridType::Codim<0>::LeafIterator Iterator;
    Iterator eendit = grid.leafend<0>();
    for (Iterator it = grid.leafbegin<0>(); it != eendit; ++it) {
        Dune::GeometryType gt = it->geometry().type();
        const Dune::FieldVector<NumberType,dim>& local = Dune::ReferenceElements<NumberType,dim>::general(gt).position(0, 0);
        Dune::FieldVector<NumberType,dim> global = it->geometry().global(local);
        if (global[0] < 1.0)
            subGridLeft.addPartial(it);
        else
            subGridRight.addPartial(it);
    }
    subGridLeft.createEnd();
    subGridRight.createEnd();

    DiffusionParameters<SubGridType,NumberType> problem;

    typedef Dune::LeafP1BoxDiffusion<SubGridType, NumberType> Diffusion;
    Diffusion diffusionLeft(subGridLeft, problem);
    Diffusion diffusionRight(subGridRight, problem);

    typedef Dune::CoupledDiffusion<Diffusion> CoupledModel;
    //typedef Dune::CoupledModel<Diffusion,Diffusion> CoupledModel;
    bool assembleGlobalMatrix = true;
    bool stationary = true;
    CoupledModel coupledModel(subGridLeft, diffusionLeft, subGridRight, diffusionRight, assembleGlobalMatrix, stationary);

    coupledModel.initial();
    coupledModel.assemble();
//    printmatrix(std::cout, coupledModel.matrix(), "global stiffness matrix", "row", 11, 4);
//    printvector(std::cout, coupledModel.rhs(), "global right hand side", "row", 200, 1, 3);
//    printvector(std::cout, coupledModel.sol(), "global solution before", "row", 200, 1, 3);
    coupledModel.solve();
//    printvector(std::cout, coupledModel.sol(), "global solution", "row", 200, 1, 3);
    coupledModel.vtkout("test_coupled", 0);

    double errorCoupled = std::max(discreteError(subGridLeft, *(*diffusionLeft), problem), discreteError(subGridRight, *(*diffusionRight), problem));
    std::cout << "error of coupled solution = " << errorCoupled << std::endl;

    DiffusionParameters<GridType,NumberType> problemGlobal;
    typedef Dune::LeafP1BoxDiffusion<GridType, NumberType> DiffusionGlobal;
    DiffusionGlobal diffusionGlobal(grid, problemGlobal);
    Dune::TimeLoop<GridType, DiffusionGlobal> timeloop(0, 1, 1, "diffusionGlobal", 1);
    timeloop.execute(diffusionGlobal);

    std::cout << "error of original global solution = " << discreteError(grid, *(*diffusionGlobal), problem) << std::endl;

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
