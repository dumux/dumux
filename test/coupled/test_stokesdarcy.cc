#include "config.h"
#include <iostream>
#define DUMMY
#ifdef DUMMY
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <../../../dune-subgrid/subgrid/subgrid.hh>
#include "dumux/stokes/dgstokes.hh"
#include "dumux/stokes/l2error.hh"
#include "dumux/stokes/h1error.hh"
#include <dumux/timedisc/timeloop.hh>
#include <dumux/coupled/coupledstokesdarcy.hh>
#include "lshapedproblem.hh"
#include "darcyproblem.hh"
#include "boxdarcy.hh"

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
    typedef MultipleCodimMultipleGeomTypeMapper<GV,P1Layout> VM;
    typedef typename GV::template Codim<dim>::Iterator VertexIterator;

    double error = 0.0;
    const GV& gridview(grid.leafView());
    VM vertexMapper(gridview);

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

int main(int argc, char** argv)
{
    try{
        const int dim=2;
        typedef double NumberType;

        // geometry
        //typedef Dune::ALUSimplexGrid<dim,dim> GridType;
        //typedef Dune::SGrid<dim,dim> GridType;
        typedef Dune::UGGrid<dim> GridType;
        Dune::GridPtr<GridType> gridPtr( argv[1] );
        GridType& grid = *gridPtr;



        // subdivide grid in subgrids
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
            if (global[0] > 1.5 && global[1] < 0.5)
                subGridRight.addPartial(it);
            else
                subGridLeft.addPartial(it);
        }
        subGridLeft.createEnd();
        subGridRight.createEnd();

        const int vOrder = 2;
        const int pOrder = 1;
        DGStokesParameters parameters;
        parameters.sigma = 0.1;
        parameters.mu = 0.01;
        Dune::LShapedProblem<SubGridType, double> stokesProblem;
        typedef Dune::DGStokes<SubGridType, vOrder, pOrder> DGStokes;
        DGStokes dGStokes(subGridLeft, stokesProblem, parameters);

        Dune::DarcyProblem<SubGridType,NumberType> darcyProblem;
        typedef Dune::LeafP1BoxDarcy<SubGridType, NumberType> DarcyModel;
        DarcyModel darcyModel(subGridRight, darcyProblem);

        typedef Dune::CoupledStokesDarcy<DGStokes,DarcyModel> CoupledModel;
        bool assembleGlobalMatrix = true;
        CoupledModel coupledModel(subGridLeft, dGStokes, subGridRight, darcyModel, assembleGlobalMatrix);

        coupledModel.initial();
        coupledModel.assemble();
        //    printmatrix(std::cout, coupledModel.matrix(), "global stiffness matrix", "row", 11, 4);
        //    printvector(std::cout, coupledModel.rhs(), "global right hand side", "row", 200, 1, 3);
        coupledModel.solve();
        //    printvector(std::cout, coupledModel.sol(), "global solution", "row", 200, 1, 3);
        //    printvector(std::cout, dGStokes.sol(), "local solution", "row", 200, 1, 3);
        coupledModel.vtkout("test_coupled", 0);


        std::cout << "Stokes L2Error velocity: ";
        for (int i = 0; i < dim; i++)
            std::cout << dGStokes.l2errorStokesSystem(i) << ", ";
        std::cout << std::endl;
        std::cout << "Stokes L2Error pressure: "<< dGStokes.l2errorStokesSystem(dim) << std::endl;
        std::cout << "Stokes H1Error velocity: ";
        for (int i = 0; i < dim; i++)
            std::cout << dGStokes.h1errorStokesSystem(i) << ", ";
        std::cout << std::endl;

        std::cout << "Darcy discrete error = " << discreteError(subGridRight, darcyModel.sol(), darcyProblem) << std::endl;

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
