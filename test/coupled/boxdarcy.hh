#ifndef DUNE_BOXDARCY_HH
#define DUNE_BOXDARCY_HH

#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include "dumux/operators/p1operatorextended.hh"
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
#include "dumux/nonlinear/nonlinearmodel.hh"
#include "dumux/fvgeometry/fvelementgeometry.hh"
#include "dumux/nonlinear/newtonmethod.hh"
#include "boxdarcyjacobian.hh"
#include "dumux/pardiso/pardiso.hh"

namespace Dune
{
template<int dim>
struct VertexLayout
{
    bool contains (Dune::GeometryType gt)
    {
        return (gt.dim() == 0);
    }
};

template<int dim>
struct ElementLayout
{
    bool contains (Dune::GeometryType gt)
    {
        return (gt.dim() == dim);
    }
};

template<class PressureFunction, class Vector, class Grid, class Problem>
void calculateDarcyVelocity(const Grid& grid, const Problem& problem, PressureFunction& pressure, Vector& xVelocity, Vector& yVelocity)
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

    ElementIterator endEIt = grid.template leafend<0>();
    for (ElementIterator eIt = grid.template leafbegin<0>(); eIt != endEIt; ++eIt)
        {
            const Element& element = *eIt;

            Dune::GeometryType geomType = element.geometry().type();

            const Dune::FieldVector<Scalar,dim>& local = Dune::ReferenceElements<Scalar,dim>::general(geomType).position(0, 0);
            Dune::FieldVector<Scalar,dim> global = element.geometry().global(local);

            int eIdx = elementMapper.map(element);

            FieldVector<Scalar, dim> pressGradient;

            for (int comp = 0; comp < dim; comp++) {
                FieldVector<int, dim> order(0);
                order[comp] = 1;
                pressGradient[comp] = pressure.derivativelocal (0, order, element, local);
            }

            FieldMatrix<Scalar, dim, dim> K = problem.K(global, element, local);

            FieldVector<Scalar, dim> darcyVelocity(0);
            K.umv(pressGradient, darcyVelocity);

            darcyVelocity *= -1.0;

            xVelocity[eIdx] = darcyVelocity[0];
            yVelocity[eIdx] = darcyVelocity[1];
        }
}

/** \todo Please doc me! */

template<class Grid, class Scalar, class ProblemType, class LocalJacobian,
         class FunctionType, class OperatorAssembler>
class BoxDarcy
    : public NonlinearModel<Grid, Scalar, ProblemType, LocalJacobian, FunctionType, OperatorAssembler>
{
public:
    typedef Dune::NonlinearModel<Grid, Scalar, ProblemType, LocalJacobian,
                                 FunctionType, OperatorAssembler> NonlinearModel;

    BoxDarcy(const Grid& grid, ProblemType& prob)
        : NonlinearModel(grid, prob), uOldTimeStep(grid)
    { }

    BoxDarcy(const Grid& grid, ProblemType& prob, int level)
        : NonlinearModel(grid, prob, level), uOldTimeStep(grid, level)
    {   }

    virtual void initial() = 0;

    virtual void update(double& dt) = 0;

    virtual void solve() = 0;

    FunctionType uOldTimeStep;
};



/** \todo Please doc me! */

template<class Grid, class Scalar, int numEq=1>
class LeafP1BoxDarcy : public BoxDarcy<Grid, Scalar, CoupledPorousMediaProblem<Grid, Scalar>, BoxDarcyJacobian<Grid, Scalar>,
                                       LeafP1Function<Grid, Scalar, numEq>, LeafP1OperatorAssembler<Grid, Scalar, numEq> >
{
public:
    // define the function type:
    typedef LeafP1Function<Grid, Scalar, numEq> FunctionType;

    // define the operator assembler type:
    typedef LeafP1OperatorAssembler<Grid, Scalar, numEq> OperatorAssembler;

    typedef CoupledPorousMediaProblem<Grid, Scalar> ProblemType;
    //    typedef DiffusionParameters<Grid, Scalar> ProblemType;

    typedef Grid GridType;

    typedef Dune::BoxDarcy<Grid, Scalar, ProblemType, BoxDarcyJacobian<Grid, Scalar>,
                           FunctionType, OperatorAssembler> BoxDarcy;

    typedef LeafP1BoxDarcy<Grid, Scalar, numEq> ThisType;

    typedef BoxDarcyJacobian<Grid, Scalar> LocalJacobian;

    // mapper: one data element per vertex
    template<int dim>
    struct P1Layout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == 0;
        }
    };

    typedef typename Grid::LeafGridView GridView;
    typedef typename GridView::IndexSet IS;
    typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,P1Layout> VertexMapper;
    typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;
    typedef typename ThisType::FunctionType::RepresentationType VectorType;
    typedef typename ThisType::OperatorAssembler::RepresentationType MatrixType;
    typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;

    LeafP1BoxDarcy (const Grid& grid, CoupledPorousMediaProblem<Grid, Scalar>& prob)
        : BoxDarcy(grid, prob), grid_(grid), vertexmapper(grid, grid.leafIndexSet()),
          size((*(this->u)).size()), uOldNewtonStep(size)
    { }

    MatrixType& matrix()
    {
        return *(this->A);
    }

    VectorType& solOldNewtonStep()
    {
        return uOldNewtonStep;
    }

    virtual void initial()
    {
        typedef typename Grid::Traits::template Codim<0>::Entity Element;
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        enum{dim = Grid::dimension};
        enum{dimworld = Grid::dimensionworld};

        const GridView& gridview(grid_.leafView());
        std::cout << "initializing solution." << std::endl;
        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
            {
                // get geometry type
                Dune::GeometryType gt = it->geometry().type();

                // get element
                const Element& element = *it;

                const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
                    sfs=Dune::LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt, 1);
                int size = sfs.size();

                for (int i = 0; i < size; i++) {
                    int globalId = vertexmapper.template map<dim>(element, sfs[i].entity());

                    // initialize cell concentration
                    (*(this->u))[globalId] = 0;
                }
            }

        // set Dirichlet boundary conditions
        for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
            {
                // get geometry type
                Dune::GeometryType gt = it->geometry().type();

                // get element
                const Element& element = *it;

                const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
                    sfs=Dune::LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt, 1);
                int size = sfs.size();

                // set type of boundary conditions
                this->localJacobian().fvGeom.update(element);
                this->localJacobian().template assembleBC<LeafTag>(element);

                IntersectionIterator endit = IntersectionIteratorGetter<Grid,LeafTag>::end(element);
                for (IntersectionIterator is = IntersectionIteratorGetter<Grid,LeafTag>::begin(element);
                     is!=endit; ++is)
                    if (is->boundary())
                        {
                            for (int i = 0; i < size; i++)
                                // handle subentities of this face
                                for (int j = 0; j < ReferenceElements<Scalar,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
                                    if (sfs[i].entity() == ReferenceElements<Scalar,dim>::general(gt).subEntity(is->numberInSelf(), 1, j, sfs[i].codim()))
                                        {
                                            if (this->localJacobian().bc(i)[0] == BoundaryConditions::dirichlet)
                                                {
                                                    // get cell center in reference element
                                                    Dune::FieldVector<Scalar,dim> local = sfs[i].position();

                                                    // get global coordinate of cell center
                                                    Dune::FieldVector<Scalar,dimworld> global = it->geometry().global(local);

                                                    int globalId = vertexmapper.template map<dim>(element, sfs[i].entity());

                                                    FieldVector<BoundaryConditions::Flags, numEq> bctype = this->problem.bctype(global, element, is, local);
                                                    if (bctype[0] == BoundaryConditions::dirichlet) {
                                                        (*(this->u))[globalId] = this->problem.g(global, element, is, local);
                                                    }
                                                    else {
                                                        std::cout << global << " is considered to be a Neumann node." << std::endl;
                                                    }
                                                }
                                        }
                        }
            }

        *(this->uOldTimeStep) = *(this->u);
        return;
    }


    virtual void update(double& dt)
    {
        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);
        NewtonMethod<Grid, ThisType> newtonMethod(grid_, *this);
        newtonMethod.execute();
        dt = this->localJacobian().getDt();
        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    virtual void solve()
    {
        Operator op(*(this->A));  // make operator out of matrix
        double red=1E-14;

#ifdef HAVE_PARDISO
        SeqPardiso<MatrixType,VectorType,VectorType> pardiso(*(this->A));
        LoopSolver<VectorType> solver(op, pardiso, red, 10, 2);
#else
        SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
        BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,1);         // an inverse operator
#endif
        InverseOperatorResult r;
        solver.apply(*(this->u), *(this->f), r);

        return;
    }

    virtual void globalDefect(FunctionType& defectGlobal) {
        typedef typename Grid::Traits::template Codim<0>::Entity Element;
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        enum {dim = Grid::dimension};
        typedef array<BoundaryConditions::Flags, numEq> BCBlockType;

        const GridView& gridview(grid_.leafView());
        (*defectGlobal)=0;

        // allocate flag vector to hold flags for essential boundary conditions
        std::vector<BCBlockType> essential(this->vertexmapper.size());
        for (typename std::vector<BCBlockType>::size_type i=0; i
                 <essential.size(); i++)
            essential[i].assign(BoundaryConditions::neumann);

        // iterate through leaf grid
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get element
            const Element& element = *it;
            this->localJacobian().fvGeom.update(element);
            int size = this->localJacobian().fvGeom.numVertices;

            this->localJacobian().setLocalSolution(element);
            this->localJacobian().computeElementData(element);
            this->localJacobian().updateVariableData(element, this->localJacobian().u);
            this->localJacobian().template localDefect<LeafTag>(element, this->localJacobian().u);

            // begin loop over vertices
            for (int i=0; i < size; i++) {
                int globalId = this->vertexmapper.template map<dim>(element,i);
                for (int equationnumber = 0; equationnumber < numEq; equationnumber++) {
                    if (this->localJacobian().bc(i)[equationnumber] == BoundaryConditions::neumann)
                        (*defectGlobal)[globalId][equationnumber]
                            += this->localJacobian().def[i][equationnumber];
                    else
                        essential[globalId].assign(BoundaryConditions::dirichlet);
                }
            }
        }

        for (typename std::vector<BCBlockType>::size_type i=0; i
                 <essential.size(); i++)
            for (int equationnumber = 0; equationnumber < numEq; equationnumber++) {
                if (essential[i][equationnumber] == BoundaryConditions::dirichlet)
                    (*defectGlobal)[i][equationnumber] = 0;
            }
    }

    virtual void vtkout (const char* name, int k)
    {
        BlockVector<FieldVector<Scalar, 1> > xVelocity(grid_.size(0));
        BlockVector<FieldVector<Scalar, 1> > yVelocity(grid_.size(0));

        calculateDarcyVelocity(grid_, this->problem, this->u, xVelocity, yVelocity);
        VTKWriter<typename Grid::LeafGridView> vtkwriter(grid_.leafView());
        vtkwriter.addVertexData(*(this->u),"pressure");
        vtkwriter.addCellData(xVelocity,"x-velocity");
        vtkwriter.addCellData(yVelocity,"y-velocity");
        vtkwriter.write(name, VTKOptions::ascii);
    }

    const Grid& grid() const
    { return grid_; }


protected:
    const Grid& grid_;
    VertexMapper vertexmapper;
    int size;
    VectorType uOldNewtonStep;
};

}
#endif
