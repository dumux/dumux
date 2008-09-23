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
#include <dune/disc/stokes/dgstokes.hh>
#include "dumux/operators/p1operatorextended.hh"
#include <dumux/timedisc/timeloop.hh>
#include <dumux/coupled/coupledstokesdarcy.hh>
#include "boxdiffusion.hh"

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;
    typedef double NumberType;

    Dune::FieldVector<double,dim> length(1.0);
    length[0] = 2.0;
    Dune::FieldVector<int,dim> size(1);
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

    int refinementSteps = 0;
    const int vOrder = 2;
    const int pOrder = 1;
    DGStokesParameters parameters;
    Example<dim, NumberType> exactSolution;
    DirichletBoundary<SubGridType> dirichletBoundary(exactSolution);
    RightHandSide<SubGridType> rightHandSide(exactSolution);
    typedef Dune::DGStokes<SubGridType, vOrder, pOrder> DGStokes;
	DGStokes dGStokes(subGridLeft, exactSolution, parameters, dirichletBoundary, rightHandSide, refinementSteps);

    DiffusionParameters<SubGridType,NumberType> problem;
    typedef Dune::LeafP1BoxDiffusion<SubGridType, NumberType> Diffusion;
    Diffusion diffusion(subGridRight, problem);

	typedef Dune::CoupledStokesDarcy<DGStokes,Diffusion> CoupledModel;
    bool assembleGlobalMatrix = true;
    CoupledModel coupledModel(subGridLeft, dGStokes, subGridRight, diffusion, assembleGlobalMatrix);

    coupledModel.initial();
    coupledModel.assemble();
    printmatrix(std::cout, coupledModel.matrix(), "global stiffness matrix", "row", 11, 4);
	printvector(std::cout, coupledModel.rhs(), "global right hand side", "row", 200, 1, 3);
	printvector(std::cout, coupledModel.sol(), "global solution before", "row", 200, 1, 3);

//    length[0] = length[1];
//    size[0] = size[1];
//    GridType controlGrid(length,size,periodic,overlap);
//    DirichletBoundary<GridType> controlDirichletBoundary(exactSolution);
//    RightHandSide<GridType> controlRightHandSide(exactSolution);
//    typedef Dune::DGStokes<GridType, vOrder, pOrder> DGStokesControl;
//	DGStokesControl dGStokesControl(controlGrid, exactSolution, parameters, controlDirichletBoundary, controlRightHandSide, refinementSteps);
//    dGStokesControl.initial();
//    dGStokesControl.assemble();
//    printmatrix(std::cout, dGStokesControl.matrix(), "Stokes stiffness matrix", "row", 11, 4);
//	printvector(std::cout, dGStokesControl.rhs(), "Stokes right hand side", "row", 200, 1, 3);
//	printvector(std::cout, dGStokesControl.sol(), "Stokes solution before", "row", 200, 1, 3);
	
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
