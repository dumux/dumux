//===============================================================
// $Id$
//===============================================================
//
// Configuration parameters

#include"config.h"                  // autoconf defines, needed by the dune headers

// select a solver
#define LOOPMGC  0
#define BICGMGC  1
#define CGBPX    0
#define CGAMG    0
#define CGILU0   0

// default parameters for adaptive algorithm
double tolerance=0;
double fraction=0.25;

// system headers
#include<iostream>               // for input/output to shell
#include<fstream>                // for input/output to files
#include<vector>                 // STL vector class
#include<sstream>
#include<complex>

#include<math.h>                 // Yes, we do some math here
#include<stdio.h>                // There is nothing better than sprintf
#include<sys/times.h>            // for timing measurements

// dune headers
#include <dune/common/fixedarray.hh>   // defines simple array classes
#include <dune/common/geometrytype.hh>
#include <dune/grid/sgrid.hh>          // a complete structured grid
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/universalmapper.hh>
#include <dune/grid/common/quadraturerules.hh>
#include <dune/common/collectivecommunication.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/vbvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>
#include <dune/istl/gsetc.hh>
#include <dune/istl/ilu.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/disc/functions/functions.hh>
#include <dune/disc/functions/p0function.hh>
#include <dune/disc/functions/p1function.hh>
#include <dune/disc/operators/p1operator.hh>
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/disc/groundwater/groundwater.hh>
#include <dune/disc/groundwater/p1groundwater.hh>
#include <dune/disc/groundwater/p1groundwaterestimator.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#if HAVE_MPI
#include <dune/istl/schwarz.hh>
#endif

#include <dune/grid/utility/gridtype.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>

// dd headers
#include"testproblem.hh"
#include"mgpreconditioner.hh"

//! Parameter for mapper class
template<int dim>
struct P0Layout
{
	bool contains (Dune::GeometryType gt)
	{
		if (gt.dim()==dim) return true;
		return false;
	} 
}; 

template<class G>
void testadaptivity (G& grid, int maxsteps, int modulo, bool globalrefine=false, bool picture=false)
{

	// first we extract the dimensions of the grid  
	const int dim = G::dimension;

	// type used for coordinates in the grid
	// such a type is exported by every grid implementation
	typedef typename G::ctype ct;
	typedef double NumberType;

	typedef Dune::LeafP1Function<G,NumberType> LeafFunction;
	typedef Dune::LeafP1OperatorAssembler<G,NumberType,1> LeafOperatorAssembler; 
	typedef Dune::GroundwaterEquationLocalStiffness<G,NumberType> LocalStiffnessType;
	typedef typename Dune::LeafP1Function<G,NumberType>::RepresentationType VectorType;
	typedef typename Dune::LeafP1OperatorAssembler<G,NumberType,1>::RepresentationType MatrixType;
	int codim = dim;

	// Get the iterator type
	// Note the use of the typename and template keywords 
	typedef typename G::template Codim<0>::LevelIterator ElementLevelIterator;

	Dune::Timer total;

	// iterate through all entities of codim 0 on the given level
	grid.globalRefine(maxsteps);
	double distance = 0.50001;
	double scaling = 0.5;

	double totaltime=0;

	// adaptation loop
	for (int step=1; step<=1; step++)
	{

		// make  functions for solution and right hand side
		LeafFunction u(grid);
		LeafFunction f(grid);
		*u = 0;

		std::cout << "=== [ STEP=" << step << " DOF=" <<  (*u).size() << std::endl;
		Dune::gridinfo(grid," ");

		typedef Dune::FieldVector<NumberType, dim> R1;
		typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<G,P0Layout> MapperType;
		typedef Dune::LeafP0Function<G,NumberType,(int)(0.5*dim*(dim+1))> KFineType;
		typedef Dune::LevelP0Function<G,NumberType,(int)(0.5*dim*(dim+1))> KC;
		MapperType mapper(grid); 
		TestModelProblem<G,NumberType> mp;

		// compute total time for each step and accumulate
		total.reset();

		LeafOperatorAssembler A(grid);

		LocalStiffnessType lstiff(mp,false);

		A.assemble(lstiff,u,f);

		// prepare solvers
		typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator; 
		Operator op(*A);  // make operator out of matrix
		double red=1E-10;
		if (step==1) red=1E-14;

		// single level preconditioner
#if CGILU0
		Dune::SeqILU0<MatrixType,VectorType,VectorType> ilu0(*A,1.0);// a precondtioner
		Dune::CGSolver<VectorType> solver(op,ilu0,red,10000,1);         // an inverse operator 
		//Dune::BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,1);         // an inverse operator 
#endif

		// geometric multigrid
#if LOOPMGC
		LocalStiffnessType mglstiff(mp,true); // with level bnd as Dirichlet
		Dune::SeqP1GeomMG<MatrixType,G,LocalStiffnessType,VectorType,VectorType,1> 
		mgc(*A,grid,mglstiff,1,0,2,false);
		Dune::LoopSolver<VectorType> solver(op,mgc,red,10000,1); 
#endif

#if BICGMGC
		LocalStiffnessType mglstiff(mp,true); // with level bnd as Dirichlet
		Dune::SeqP1GeomMG<MatrixType,G,LocalStiffnessType,VectorType,VectorType,1> 
		mgc(*A,grid,mglstiff,0,1,1,false);
		Dune::CGSolver<VectorType> solver(op,mgc,red,10000,2); 
#endif

#if CGBPX
		LocalStiffnessType mglstiff(mp,true); // with level bnd as Dirichlet
		Dune::SeqP1GeomMG<MatrixType,G,LocalStiffnessType,VectorType,VectorType,1> 
		bpx(*A,grid,mglstiff,1,1,0,true);
		Dune::CGSolver<VectorType> solver(op,bpx,red,10000,1);         // an inverse operator 
#endif

		// algebraic multigrid
#if CGAMG
		typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixType,
		Dune::Amg::FirstDiagonal> > Criterion;
		typedef Dune::SeqSSOR<MatrixType,VectorType,VectorType> Smoother;
		typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
		SmootherArgs smootherArgs;
		smootherArgs.iterations = 2;
		int maxlevel = 20, coarsenTarget = 100;
		Criterion criterion(maxlevel, coarsenTarget);
		criterion.setMaxDistance(2);
		typedef Dune::Amg::AMG<Operator,VectorType,Smoother> AMG;
		AMG amg(op,criterion,smootherArgs,1,1);								       
		Dune::CGSolver<VectorType> solver(op,amg,red,10000,1);
#endif

		// solve the linear system
		Dune::InverseOperatorResult r;
		//*u = 0;
		solver.apply(*u,*f,r);		

		std::cout << "=== TIME for assembly and solve " << total.elapsed() << " second(s)" << std::endl; 

		// 	  double norm1 = L2Error(grid, u, LinearFunction<G,double,dim>(), 5);
		//  	  double norm2 = H1Error(grid, u, LinearFunction<G,double,dim>(), 5);
		// 	  std::cout << std::setw(8) << grid.size(codim) << "   L2 "
		// 		    << std::scientific << std::showpoint << std::setprecision(10) << " " << norm1;
		// 	  std::cout << std::setw(8) << "H1 " << std::scientific << std::showpoint << std::setprecision(10)
		// 		    << " " << norm2 << "  time  " << watch.elapsed() << std::endl;

		if (picture)
		{
			// graphics output
			std::ostringstream os;
			os << "u." << grid.name() << "." << dim << "d";
			if (globalrefine)
				os << ".global";
			else
				os << ".local";
			Dune::LeafVTKWriter<typename G::LeafGridView> vtkwriter(grid.leafView());
			vtkwriter.addVertexData(*u,"solution");
			std::string s(os.str());
			vtkwriter.write(s.c_str(),Dune::VTKOptions::ascii);
		}

	}
}


// the main program
int main (int argc , char ** argv)
{
	try {
		if (argc!=2)
		{
			std::cout << "usage: test_multigrid #steps" << std::endl;
			return 1;
		}
		std::string arg2(argv[1]);

		std::istringstream is2(arg2);
		int steps;
		is2 >> steps;

			bool globalrefine=true;

			bool picture=false;

			const int dim = 2;

			typedef Dune::SGrid<dim,dim> GridType;
			//typedef Dune::UGGrid<dim> GridType;

			// use unitcube from grids 
			std::stringstream dgfFileName;
			dgfFileName << "grids/unitcube" << GridType :: dimension << ".dgf";

			// create grid pointer, GridType is defined by gridtype.hh
			Dune::GridPtr<GridType> gridPtr( dgfFileName.str() );

			// grid reference 
			GridType& grid = *gridPtr;

			testadaptivity(grid, steps, 1, globalrefine, picture);
	}
	catch (Dune::Exception& error)
	{
		std::cout << error << std::endl;
	}

	return 0;

}
