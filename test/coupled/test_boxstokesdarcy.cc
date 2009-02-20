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

template<class Vector, class Grid, class Problem>
void discreteStokesError(const Grid& grid, const Problem& problem, Vector& solution)
{
    typedef typename Grid::ctype Scalar;
    enum {dim=Grid::dimension};
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Grid::LeafGridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename Grid::LeafGridView::IndexSet IS;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IS,ElementLayout> ElementMapper;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IS,VertexLayout> VertexMapper;

    VertexMapper vertexMapper(grid, grid.leafView().indexSet());
    ElementMapper elementMapper(grid, grid.leafView().indexSet());

    Scalar errPressure = 0;
    Scalar errVelocity = 0;
    Scalar constant = 0;
    ElementIterator endEIt = grid.template leafend<0>();
    for (ElementIterator eIt = grid.template leafbegin<0>(); eIt != endEIt; ++eIt)
    {
        const Element& element = *eIt;

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
            if (global[0] < 1)
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

        typedef Dune::BoxStokesDarcy<BoxStokes,DarcyModel,NumberType> CoupledModel;
        bool assembleGlobalMatrix = true;
        CoupledModel coupledModel(subGridStokes, boxStokes, subGridDarcy, darcyModel, assembleGlobalMatrix);

        Dune::TimeLoop<GridType, CoupledModel, false> timeloop(0, 1, 1, "test_boxstokesdarcy", 1);
        Dune::Timer timer;
        timer.reset();
        timeloop.execute(coupledModel);
        std::cout << "timeloop.execute took " << timer.elapsed() << " seconds" << std::endl;

//        printvector(std::cout, coupledModel.sol(), "global solution", "row", 200, 1, 3);
//        printvector(std::cout, boxStokes.sol(), "local Stokes solution", "row", 200, 1, 3);
//        printvector(std::cout, darcyModel.sol(), "local Darcy solution", "row", 200, 1, 3);

        std::cout << "Darcy discrete error = " << discreteDarcyError(subGridDarcy, darcyModel.sol(), darcyParam) << std::endl;
        discreteStokesError(subGridStokes, stokesProblem, boxStokes.u);

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
