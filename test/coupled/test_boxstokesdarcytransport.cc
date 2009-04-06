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
#include <dumux/2p2c/fv/box2p2c.hh>
#include "lshapedproblem.hh"
#include "../stokes/boxstokes.hh"
#include <dumux/coupled/boxstokesdarcytransport.hh>

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

template<class Vector, class Grid, class Problem>
void discreteStokesError(const Grid& grid, const Problem& problem, Vector& solution)
{
    typedef typename Grid::ctype Scalar;
    enum {dim=Grid::dimension};
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Grid::LeafGridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename Grid::LeafGridView::IndexSet IS;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<typename Grid::LeafGridView,ElementLayout> ElementMapper;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<typename Grid::LeafGridView,VertexLayout> VertexMapper;

    VertexMapper vertexMapper(grid.leafView());
    ElementMapper elementMapper(grid.leafView());

    Scalar errPressure = 0;
    Scalar errVelocity = 0;
    Scalar constant = 0;
    ElementIterator endEIt = grid.template leafend<0>();
    for (ElementIterator elementIt = grid.template leafbegin<0>(); elementIt != endEIt; ++elementIt)
    {
        const Element& element = *elementIt;

        Dune::GeometryType geomType = element.geometry().type();

        const Dune::FieldVector<Scalar,dim>& local = Dune::ReferenceElements<Scalar,dim>::general(geomType).position(0, 0);
        Dune::FieldVector<Scalar,dim> global = element.geometry().global(local);

        Scalar volume = element.geometry().integrationElement(local)
            *Dune::ReferenceElements<Scalar,dim>::general(geomType).volume();

        int eIdx = elementMapper.map(element);

        Scalar approxPressure = solution.evallocal (2, element, local);
        Scalar exactPressure = problem.pressure(global);
        if (eIdx == 0)
            constant = exactPressure - approxPressure;
        approxPressure += constant;
        errPressure += volume*(approxPressure - exactPressure)*(approxPressure - exactPressure);

        Scalar approxXV = solution.evallocal (0, element, local);
        Scalar exactXV = problem.velocity(global)[0];
        errVelocity += volume*(approxXV - exactXV)*(approxXV - exactXV);

        Scalar approxYV = solution.evallocal (1, element, local);
        Scalar exactYV = problem.velocity(global)[1];
        errVelocity += volume*(approxYV - exactYV)*(approxYV - exactYV);
    }

    errPressure = sqrt(errPressure);
    errVelocity = sqrt(errVelocity);

    std::cout << "Error in discrete L2 norm:\nPressure: " << errPressure << "\nVelocity: " << errVelocity << std::endl;
}

template<class Grid, class Solution, class Problem>
double discreteDarcyError(const Grid& grid, const Solution& solution, const Problem& problem)
{
    enum{dim=Grid::dimension};
    typedef typename Grid::LeafGridView GV;
    typedef typename GV::IndexSet IS;
    typedef MultipleCodimMultipleGeomTypeMapper<GV,P1Layout> VM;
    typedef typename GV::template Codim<dim>::Iterator VertexIterator;

    double error = 0.0;
    const GV& gridview(grid.leafView());
    VM vertexMapper(gridview);

    VertexIterator endIt = gridview.template end<dim>();
    VertexIterator vertexIt = gridview.template begin<dim>();
    for (; vertexIt != endIt; ++vertexIt)
    {
        // get exact solution at vertex
        FieldVector<double,dim> globalCoord = (*vertexIt).geometry().corner(0);
        double exact = problem.exact(globalCoord);

        // get approximate solution at vertex
        int globalId = vertexMapper.map(*vertexIt);
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
        typedef double Scalar;

        // grid and geometry
        //typedef Dune::ALUSimplexGrid<dim,dim> GridType;
        typedef Dune::SGrid<dim,dim> GridType;

        if (argc != 2) {
            std::cout << boost::format("usage: %s grid\n")%argv[0];
            return 1;
        }
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
        for (Iterator elementIt = grid.leafbegin<0>(); elementIt != eendit; ++elementIt) {
            Dune::GeometryType gt = elementIt->geometry().type();
            const Dune::FieldVector<Scalar,dim>& local = Dune::ReferenceElements<Scalar,dim>::general(gt).position(0, 0);
            Dune::FieldVector<Scalar,dim> global = elementIt->geometry().global(local);
            if (global[0] < 1.5 || global[1] > 0.5)
                subGridStokes.insert(*elementIt);
            else
                subGridDarcy.insert(*elementIt);
        }
        subGridStokes.createEnd();
        subGridDarcy.createEnd();

        // choose fluids and properties
        Dune::Liq_WaterAir wPhase;
        Dune::Gas_WaterAir nPhase;
        Dune::TwoPTwoCDarcySoil<SubGridType, Scalar> soil;

        Dune::TwoPhaseRelations<SubGridType, Scalar> materialLaw(soil, wPhase, nPhase);
        Dune::CWaterAir multicomp(wPhase, nPhase);
        Scalar depthBOR = 1.0;

        Dune::LShapedProblem<SubGridType, Scalar> stokesProblem;
        typedef Dune::LeafP1BoxStokes<SubGridType, Scalar, dim> StokesModel;
        StokesModel stokesModel(subGridStokes, stokesProblem);

        Dune::TwoPTwoCDarcyProblem<SubGridType,Scalar> darcyProblem(wPhase, nPhase, soil,
														depthBOR, materialLaw, multicomp);
        // create two-phase two-component problem
        typedef Dune::VtkMultiWriter<GridType::LeafGridView> MultiWriter;
        typedef Dune::Box2P2C<SubGridType, Scalar, MultiWriter> DarcyModel;
//        typedef Dune::LeafP1BoxDarcyTransport<SubGridType, Scalar> DarcyModel;
        DarcyModel darcyModel(subGridDarcy, darcyProblem);

        typedef Dune::BoxStokesDarcyTransport<StokesModel, DarcyModel, Scalar> CoupledModel;
        bool assembleGlobalMatrix = true;
        CoupledModel coupledModel(subGridStokes, stokesModel, subGridDarcy, darcyModel, assembleGlobalMatrix);

        coupledModel.vtkout("initial", 0);

        Dune::TimeLoop<GridType, CoupledModel, false> timeloop(0, 1, 1, "test_boxstokesdarcy", 1);
        Dune::Timer timer;
        timer.reset();
        timeloop.execute(coupledModel);
        std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

        //        printvector(std::cout, coupledModel.sol(), "global solution", "row", 200, 1, 3);
        //        printvector(std::cout, boxModel.sol(), "local Stokes solution", "row", 200, 1, 3);
        //        printvector(std::cout, darcyModel.sol(), "local Darcy solution", "row", 200, 1, 3);

//        std::cout << "Darcy discrete error = " << discreteDarcyError(subGridDarcy, darcyModel.sol(), darcyProblem) << std::endl;
        discreteStokesError(subGridStokes, stokesProblem, stokesModel.u);

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
