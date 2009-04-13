#include "config.h"
#include <iostream>
#include <boost/format.hpp>
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
#include <dumux/timedisc/timeloop.hh>
#include "dumux/io/vtkmultiwriter.hh"

#include "2p2cdarcyproblem.hh"
#include "boxdarcytransport.hh"
#include "2cstokesproblem.hh"
#include "boxstokestransport.hh"
#include <dumux/coupled/boxstokesdarcytransport.hh>

//namespace Dune
//{
//template<int dim>
//struct P1Layout
//{
//    bool contains (Dune::GeometryType gt)
//    {
//        return gt.dim() == 0;
//    }
//};

//template<int dim>
//struct NodeLayout
//{
//    bool contains(Dune::GeometryType gt) {
//        return gt.dim() == 0;
//    }
//};

int main(int argc, char** argv)
{
    try{
        const int dim=2;
        typedef double Scalar;

        // grid and geometry
        //typedef Dune::ALUSimplexGrid<dim,dim> Grid;
        typedef Dune::SGrid<dim,dim> Grid;
        typedef Grid::LeafGridView GridView;
        typedef Dune::VtkMultiWriter<GridView> MultiWriter;

        if (argc != 2) {
            std::cout << boost::format("usage: %s grid\n")%argv[0];
            return 1;
        }
        Dune::GridPtr<Grid> gridPtr( argv[1] );
        Grid& grid = *gridPtr;

        // subdivide grid in subgrids
        typedef Dune::SubGrid<dim,Grid> SubGrid;
        SubGrid subGridStokes(grid);
        SubGrid subGridDarcy(grid);
        subGridStokes.createBegin();
        subGridDarcy.createBegin();
        typedef Grid::Codim<0>::LeafIterator Iterator;

        Iterator eendit = grid.leafend<0>();
        for (Iterator elementIt = grid.leafbegin<0>(); elementIt != eendit; ++elementIt) {
            Dune::GeometryType gt = elementIt->geometry().type();
            const Dune::FieldVector<Scalar,dim>& local = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);
            Dune::FieldVector<Scalar,dim> global = elementIt->geometry().global(local);
            if (global[0] < 1.5 || global[1] > 0.5 || global[0] > 4.0)
                subGridStokes.insert(*elementIt);
            else
                subGridDarcy.insert(*elementIt);
        }
        subGridStokes.createEnd();
        subGridDarcy.createEnd();

        // choose fluids and properties
        Dune::Liq_WaterAir wPhase;
        Dune::Gas_WaterAir nPhase;
        Dune::TwoPTwoCDarcySoil<SubGrid, Scalar> soil;

        // create twophase and multicomponent relations
        Dune::TwoPhaseRelations<SubGrid, Scalar> materialLaw(soil, wPhase, nPhase);
        Dune::CWaterAir multicomp(wPhase, nPhase);
        Scalar depthBOR = 1.0;

        Dune::TwoCStokesProblem<SubGrid, Scalar> stokesProblem(nPhase, soil, multicomp);
        typedef Dune::LeafP1BoxStokesTransport<SubGrid, Scalar, dim> StokesModel;
        StokesModel stokesModel(subGridStokes, stokesProblem);

        // create two-phase two-component problem and (dummy) multiwriter
        Dune::TwoPTwoCDarcyProblem<SubGrid,Scalar> darcyProblem(wPhase, nPhase, soil,
								depthBOR, materialLaw, multicomp);
        typedef Dune::BoxDarcyTransport<SubGrid, Scalar, MultiWriter> DarcyModel;
        DarcyModel darcyModel(subGridDarcy, darcyProblem);

        typedef Dune::BoxStokesDarcyTransport<StokesModel, DarcyModel, Scalar> CoupledModel;
        bool assembleGlobalMatrix = true;
        CoupledModel coupledModel(subGridStokes, stokesModel, subGridDarcy, darcyModel, assembleGlobalMatrix);
//        coupledModel.vtkout("initial", 0);

        Dune::TimeLoop<Grid, CoupledModel, false> timeloop(0, 1, 1, "test_boxcoupledtransport", 1);
        Dune::Timer timer;
        timer.reset();
        timeloop.execute(coupledModel);
        std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

        //        printvector(std::cout, coupledModel.sol(), "global solution", "row", 200, 1, 3);
        //        printvector(std::cout, boxModel.sol(), "local Stokes solution", "row", 200, 1, 3);
        //        printvector(std::cout, darcyModel.sol(), "local Darcy solution", "row", 200, 1, 3);

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
