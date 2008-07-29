#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include"dune/common/mpihelper.hh" // An initializer of MPI
#include"dune/common/exceptions.hh" // We use exceptions
#include "dune_multiscale_headers.hh"
#include <dune/grid/uggrid.hh>
#include <dune/grid/albertagrid.hh>
#include <dune/grid/alugrid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/grid/sgrid.hh> // load sgrid definition

#include <dumux/timedisc/timeloop.hh>
//#include <dumux/material/properties.hh>
#include "boxdiffusion_couple.hh"
// timeloop options
void inline TimeloopOptsPipe( double& tstart, double& tend, double& max_dt, double& first_dt, double& CFL_factor, int& flag,
		int& n_iter, double& max_def, int& modulo, int& stages )
{
	tstart      = 0.0;      // start time of simulation
	tend        = 0.1;    // end time of simulation
	max_dt      = 0.1;      // maximum time step
	first_dt    = 1;      // maximum first time step
	CFL_factor  = 1e0;     // safety-factor for Courant-Friedrich-Levy criterion
	flag        = 0;        // 0: no iteration
	// 1: iteration with n_iter steps
	// 2: iteration with max_def defect and timestep-refinement after n_iter iterations
	n_iter      = 50;
	max_def     = 1e-12;    // tolerance
	modulo      = 1;        // vtk - output is written for every modulo-th timestep

	stages      = 0;        // = 0: implicit Euler
	// > 0: number of stages for explicit RK method
}

int main(int argc, char** argv)
{
	try{
		//Maybe initialize Mpi
		Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
		std::cout << "Hello World! This is pipe_box." << std::endl;
		if(Dune::MPIHelper::isFake)
			std::cout<< "This is a sequential program." << std::endl;
		else
			std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()
			<<" processes!"<<std::endl;

		// set up the grid
		const int dim = 3;
		typedef double NumberType;

		//      typedef Dune::ALUSimplexGrid<dim,dim> GridType;
		//      typedef Dune::UGGrid<dim> GridType;
		typedef Dune::ALUCubeGrid<dim,dim> GridType;

		// create grid pointer, GridType is defined by gridtype.hh
		Dune::GridPtr<GridType> gridPtr( "grids/csp_debug_big.dgf" );
		//    Dune::GridPtr<GridType> gridPtr( "grids/csp_dx_40cm.dgf" );

		// grid reference
		GridType& grid = *gridPtr;

		typedef VertexOnLine<GridType> VertexType1;
		typedef VertexOutLine<GridType> VertexType2;

		typedef std::vector<VertexType1> VertexVectorOnLineType;
		VertexVectorOnLineType vertexVectorOnLine;
		typedef std::vector<VertexType2> VertexVectorOutLineType;
		VertexVectorOutLineType vertexVectorOutLine;

		typedef GridType::Traits::LeafIndexSet IS;
		typedef Dune::MultipleCodimMultipleGeomTypeMapper<GridType, IS, PDimLayout > VM; //vertexmapper maps (elementpointer,localIdofnode) to globalnodeindexset
		VM vertexmapper(grid, grid.leafIndexSet());

		typedef std::map<unsigned, unsigned, GlobalNodeIdCompare> MapGlobalNodeIDtoPipeNodeOnOutlineIndexType;
		MapGlobalNodeIDtoPipeNodeOnOutlineIndexType mapGlobalNodeIDtoPipeNodeOnlineIndex;
		MapGlobalNodeIDtoPipeNodeOnOutlineIndexType mapGlobalNodeIDtoPipeNodeOutlineIndex;

		isolate(grid, gridPtr, vertexmapper, mapGlobalNodeIDtoPipeNodeOnlineIndex, mapGlobalNodeIDtoPipeNodeOutlineIndex, vertexVectorOnLine, vertexVectorOutLine);

		typedef Dune::FieldVector<double,GridType::dimension>  FieldVector;

		typedef PressureBoundary BCP;
		typedef VelocityBoundary BCV;

		// set the initial condition
		typedef ICPressurePipe<FieldVector> ICP;
		typedef ICVelocityPipe<FieldVector> ICV;

		// set the initial condition
		typedef sourceSinkPipe<FieldVector> SST;

		// set the mobility/saturation relation
		//typedef LinearPressure<Air, Air> Mob ;

		// set fluid property
		Water WaterProp;

		// set the numerical flux function
		typedef Lambda<Water> Lmbd;

		// set the block vector type for pressure
		typedef Dune::BlockVector< Dune::FieldVector<double,1> > PresType;

		// make a mapper for codim 0 entities in the leaf grid
		typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,P0Layout> MapperTypeC0;
		MapperTypeC0 mapperC0(grid);

		// make a mapper for codim DIM entities in the leaf grid
		typedef Dune::LeafMultipleCodimMultipleGeomTypeMapper<GridType,PDimLayout> MapperTypeCDim;
		MapperTypeCDim mapperCDim(grid);

		int numNodesPorousMesh = mapperCDim.size();

		std::cout<<"number of numNodesPorousMesh: "<<numNodesPorousMesh<<std::endl;

		// calculate the number of elements
		int nElem = vertexVectorOnLine.size();
		std::cout<<"number of cells: "<<nElem<<std::endl;

		double temp = 283.15;
		double density = WaterProp.density(temp, 1.0e5);
		double kinematicViscosity = WaterProp.viscosity()/density;
		double roughness = 0.3;
		double diameter = 0.4; // in [m]


		double alphaEX;
		PresType pressurePorous(numNodesPorousMesh);
		PresType pressurePipe;


		alphaEX = 0;
		pressurePorous = 0;
		double toler = 1.0e-10;
		double diff=1e100;
		PresType oldPipePress(nElem);
		PresType oldPorousPress(numNodesPorousMesh);
		oldPipePress = 0;
		oldPorousPress = 0;
		int maxIt = 10;
		int iter = 0;
		while (diff >toler)
		{
			iter++;
			if (iter > maxIt) {
				std::cout << "Did not converge!" << std::endl;
				break;
			}

			PipeFlow< BCP, BCV, ICP, ICV, SST, PresType, Lmbd, GridType, MapperTypeC0, MapGlobalNodeIDtoPipeNodeOnOutlineIndexType, VertexVectorOnLineType, VertexVectorOutLineType> 
			pipeflow(nElem, &grid, &mapperC0, &mapGlobalNodeIDtoPipeNodeOnlineIndex, &mapGlobalNodeIDtoPipeNodeOutlineIndex, &vertexVectorOnLine, &vertexVectorOutLine,
					pressurePorous, alphaEX, temp, density, kinematicViscosity, roughness, diameter );

			pipeflow.timeLoop();

			alphaEX = 1.0e-3/(WaterProp.density() * 9.81);
			pressurePipe = pipeflow.pressure;

			typedef DiffusionParameters<GridType, NumberType, MapGlobalNodeIDtoPipeNodeOnOutlineIndexType, VM, VertexVectorOnLineType> DiffusionParameters; 
			DiffusionParameters diffusionProblem(alphaEX, pressurePipe, &vertexVectorOnLine,
					mapGlobalNodeIDtoPipeNodeOnlineIndex, vertexmapper, WaterProp);

			typedef Dune::LeafP1BoxDiffusion<GridType, NumberType, MapGlobalNodeIDtoPipeNodeOnOutlineIndexType, VM, VertexVectorOnLineType> Diffusion;
			Diffusion diffusion(grid, diffusionProblem);

			Dune::TimeLoop<GridType, Diffusion> timeloop(0, 1, 1, "box3d", 1);

			timeloop.execute(diffusion);

			pressurePorous = *(*diffusion);
			pipeflow.alphaEX = alphaEX;

			PresType difVec1 = pipeflow.pressure;
			difVec1 -= oldPipePress;
			double dif1 = difVec1.two_norm()/pipeflow.pressure.two_norm();
			PresType difVec2 = *(*diffusion);
			difVec2 -= oldPorousPress;
			double dif2 = difVec2.two_norm()/(*(*diffusion)).two_norm();
			diff = std::max(dif1, dif2);
			std::cout << "diff1 = " << dif1 << ", dif2 = " << dif2 << std::endl;

			oldPipePress = pipeflow.pressure; 
			oldPorousPress = *(*diffusion);

		}

		printvector(std::cout,oldPipePress,"PressurePipe","row",100,1,4);
//		printvector(std::cout,pressurePorous,"PressurePorous","row",100,1,4);


		//    pipeflow.PrintVertexVector();

		return 0;
	}
	catch (Dune::Exception &e){
		std::cerr << "Dune reported error: " << e << std::endl;
	}
	catch (...){
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}

