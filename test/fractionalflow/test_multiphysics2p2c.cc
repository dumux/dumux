#define DUNE_DEVEL_MODE
#include "config.h"
#include <iostream>
#include <iomanip>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/io.hh>
#include <dune/common/timer.hh>
#include "dumux/mystuff/timedisc/timeloop.hh"
#include "dumux/mystuff/material/phaseproperties/phaseproperties_waterair.hh"
#include <dumux/mystuff/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>
#include "dumux/mystuff//transport/problems/testproblem_2p2c.hh"
#include "multiphysics2p2c.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include "dumux/timedisc/expliciteulerstep.hh"

//  @author Jochen Fritz
//	last change: 19.11.08
//

int main(int argc, char** argv)
{
  try{
    // define the problem dimensions
    const int dim=2;

    // create a grid object
    typedef double NumberType;
    typedef Dune::SGrid<dim,dim> GridType;
    typedef Dune::FieldVector<GridType::ctype,dim> FieldVector;
    Dune::FieldVector<int,dim> N(100); N[0] = 100;
    FieldVector L(0);
    FieldVector H(300); H[0] = 300;
    GridType grid(N,L,H);

    grid.globalRefine(0);

    typedef GridType::Codim<0>::LevelIterator HostIterator;
    typedef GridType::Codim<0>::HierarchicIterator HierarchicIterator;
    typedef Dune::IntersectionIteratorGetter<GridType,Dune::LevelTag>::IntersectionIterator HostIntersectionIterator;

    Dune::SubGrid<dim,GridType> subGrid(grid);

		subGrid.createBegin();
			HostIterator endit = grid.lend<0>(0);
			for (HostIterator it = grid.lbegin<0>(0); it!=endit; ++it)
			{
				Dune::GeometryType gt = it->geometry().type();
				const Dune::FieldVector<GridType::ctype,dim>&
					local = Dune::ReferenceElements<GridType::ctype,dim>::general(gt).position(0,0);
				Dune::FieldVector<GridType::ctype,dim> global = it->geometry().global(local);

				if (global[0] > 45 && global[0] < 180 && global[1] > 105 && global[1] <195)
				{
					subGrid.addPartial(it);
					HierarchicIterator hend = it->hend(0);
					for ( HierarchicIterator h = it->hbegin(0); h != hend; ++h)
						subGrid.add(h);
				}
			}
		subGrid.createEnd();

    double tStart = 0;
    double tEnd = 5e6;
    int modulo = 1;
    double cFLFactor = 0.9;

    Dune::Liq_WaterAir wetmat;
    Dune::Gas_WaterAir nonwetmat;
//    Dune::HomogeneousSoil<GridType, NumberType> soil;
    Dune::HeterogeneousSoil<GridType, NumberType> soil(grid, "permeab.dat", false);
    soil.permeability.vtkout("perm",grid);

//    std::ofstream outf("permeab.dat");
//    for(int i = 0;i<grid.size(0); i++)
//    {
//    	outf<<(*(soil.permeability))[i]<<","<<std::endl;
//    }

    Dune::TwoPhaseRelations<GridType, NumberType> materialLaw(soil, wetmat, nonwetmat);

    Dune::VariableClass2p2c<GridType,NumberType> var(grid);

    typedef Dune::Testproblem_2p2c<GridType, NumberType> TransProb;
    TransProb problem(grid, var, wetmat, nonwetmat, soil, grid.maxLevel(), materialLaw, false);

    Dune::DiffusivePart<GridType, NumberType> diffPart;
    const Dune::Upwind<NumberType> numFl;

    typedef Dune::Multiphysics2p2c<GridType, NumberType> ModelType;
    ModelType model(grid, subGrid, problem, grid.maxLevel(), diffPart, false, 0.8, numFl);

    Dune::ExplicitEulerStep<GridType, ModelType> timestep;
    Dune::TimeLoop<GridType, ModelType > timeloop(tStart, tEnd, "mp", modulo, cFLFactor, 1e100, 1e100, timestep);

    Dune::Timer timer;
    timer.reset();
    timeloop.execute(model);
    std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
