#include "config.h"
#include <iostream>
#include <iomanip>

#include <dune/grid/sgrid.hh>
#include <dune/istl/io.hh>
#include <dune/istl/solvers.hh>
#include "dumux/stokes/localmixed.hh"
#include "dumux/operators/mixedoperator.hh"
#include "dumux/pardiso/pardiso.hh"

#include "yxproblem.hh"

template<int dim>
struct ElementAndFaceLayout
{
    bool contains (Dune::GeometryType gt)
    {
        return (gt.dim() == dim || gt.dim() == dim-1);
    }
};

template<class Vector, class Grid>
void calculatePressureAndVelocities(const Vector& u, const Grid& grid, Vector& pressure, Vector& xV, Vector& yV)
{

    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Grid::LeafGridView::IndexSet IS;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IS,ElementAndFaceLayout> ElementAndFaceMapper;

    ElementAndFaceMapper elementAndFaceMapper(grid, grid.leafView().indexSet());

    ElementIterator endEIt = grid.template leafend<0>();
    for (ElementIterator eIt = grid.template leafbegin<0>(); eIt != endEIt; ++eIt)
        {
            const Element& element = *eIt;

            int eIdx = elementAndFaceMapper.map(element);
            int leftFaceIdx = elementAndFaceMapper.template map<1>(element, 0);
            int rightFaceIdx = elementAndFaceMapper.template map<1>(element, 1);
            int bottomFaceIdx = elementAndFaceMapper.template map<1>(element, 2);
            int topFaceIdx = elementAndFaceMapper.template map<1>(element, 3);

            pressure[eIdx] = u[eIdx];
            xV[eIdx] = 0.5*(u[leftFaceIdx] + u[rightFaceIdx]);
            yV[eIdx] = 0.5*(u[bottomFaceIdx] + u[topFaceIdx]);
        }
}

template<class Vector, class Grid, class Problem>
void calculateError(const Grid& grid, const Problem& problem, Vector& pressure, Vector& xV, Vector& yV)
{
    typedef typename Grid::ctype Scalar;
    enum {dim=Grid::dimension};
    typedef typename Grid::template Codim<0>::Entity Element;
    typedef typename Grid::LeafGridView::template Codim<0>::Iterator ElementIterator;
    typedef typename Grid::LeafGridView::IndexSet IS;
    typedef Dune::MultipleCodimMultipleGeomTypeMapper<Grid,IS,ElementAndFaceLayout> ElementAndFaceMapper;

    ElementAndFaceMapper elementAndFaceMapper(grid, grid.leafView().indexSet());

    Scalar errPressure = 0;
    Scalar errVelocity = 0;
    Scalar constant;
    ElementIterator endEIt = grid.template leafend<0>();
    for (ElementIterator eIt = grid.template leafbegin<0>(); eIt != endEIt; ++eIt)
        {
            const Element& element = *eIt;

            Dune::GeometryType geomType = element.geometry().type();

            const Dune::FieldVector<Scalar,dim>& local = Dune::ReferenceElements<Scalar,dim>::general(geomType).position(0, 0);
            Dune::FieldVector<Scalar,dim> global = element.geometry().global(local);

            Scalar volume = element.geometry().integrationElement(local)
                *Dune::ReferenceElements<Scalar,dim>::general(geomType).volume();

            int eIdx = elementAndFaceMapper.map(element);

            if (eIdx == 0)
                constant = problem.pressure(global) - pressure[eIdx];

            Scalar approxPressure = pressure[eIdx] + constant;
            Scalar exactPressure = problem.pressure(global);
            errPressure += volume*(approxPressure - exactPressure)*(approxPressure - exactPressure);

            Scalar approxXV = xV[eIdx];
            Scalar exactXV = problem.velocity(global)[0];
            errVelocity += volume*(approxXV - exactXV)*(approxXV - exactXV);

            Scalar approxYV = yV[eIdx];
            Scalar exactYV = problem.velocity(global)[1];
            errVelocity += volume*(approxYV - exactYV)*(approxYV - exactYV);
        }

    errPressure = sqrt(errPressure);
    errVelocity = sqrt(errVelocity);

    std::cout << "Error in discrete L2 norm:\nPressure: " << errPressure << "\nVelocity: " << errVelocity << std::endl;
}

int main(int argc, char** argv)
{
    try
        {
            const int dim = 2;

            typedef double Scalar;
            typedef Dune::SGrid<dim, dim> Grid;

            typedef Dune::FieldVector<Scalar,dim> FieldVector;
            Dune::FieldVector<int,dim> N(16); N[0] = 16;
            FieldVector L(0);
            FieldVector H(1);
            Grid grid(N,L,H);

            typedef Dune::LeafMixedFunction <Grid, Scalar, 1> MixedFunction;
            MixedFunction u(grid);
            MixedFunction f(grid);

            typedef Dune::YXProblem<Grid, Scalar> Problem;
            Problem problem;

            typedef Dune::LocalMixed<Grid, Scalar, 1> LocalMixed;
            LocalMixed localMixed(problem);

            typedef Dune::LeafMixedOperatorAssembler<Grid, Scalar, 1> MixedAssembler;
            MixedAssembler mixedAssembler(grid);
            mixedAssembler.assemble(localMixed, u, f);

            //      printmatrix(std::cout, *mixedAssembler, "global stiffness matrix", "row", 11, 4);
            //      printvector(std::cout, *f, "right hand side", "row", 200, 1, 3);

            typedef Dune::BlockVector<Dune::FieldVector<Scalar, 1> > Vector;
            typedef Dune::BCRSMatrix<Dune::FieldMatrix<Scalar, 1, 1> > Matrix;
            Dune::MatrixAdapter<Matrix,Vector,Vector> op(*mixedAssembler);
            Dune::InverseOperatorResult r;
            Dune::SeqPardiso<Matrix,Vector,Vector> preconditioner(*mixedAssembler);
            Dune::LoopSolver<Vector> solver(op, preconditioner, 1E-14, 1000, 1);
            solver.apply(*u, *f, r);

            //      printvector(std::cout, *u, "solution", "row", 200, 1, 3);

            Vector pressure(grid.size(0));
            Vector xCellVelocity(grid.size(0));
            Vector yCellVelocity(grid.size(0));
            calculatePressureAndVelocities(*u, grid, pressure, xCellVelocity, yCellVelocity);

            calculateError(grid, problem, pressure, xCellVelocity, yCellVelocity);

            Dune::VTKWriter<Grid::LeafGridView> vtkwriter(grid.leafView());
            vtkwriter.addCellData(pressure, "pressure");
            vtkwriter.addCellData(xCellVelocity, "velocity x-comp");
            vtkwriter.addCellData(yCellVelocity, "velocity y-comp");
            vtkwriter.write("test_mixed", Dune::VTKOptions::ascii);

            //      typedef Dune::StokesJacobian<Grid, Scalar> StokesJacobian;
            //      typedef Dune::StokesProblem<Grid, Scalar> StokesProblem;
            //      typedef Dune::MixedNonlinearModel<Grid, Scalar, StokesProblem, StokesJacobian> NonlinearModel;
            //      NonlinearModel nonlinearModel(grid, problem);

            return 0;
        }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
