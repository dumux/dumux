#ifndef DUNE_BOXSTOKES_HH
#define DUNE_BOXSTOKES_HH

#include <dune/disc/shapefunctions/lagrangeshapefunctions.hh>
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
#include "boxstokesjacobian.hh"
#include "dumux/stokes/stokesproblem.hh"
#include "dumux/pardiso/pardiso.hh"

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar, class ProblemType, class LocalJacobian,
         class FunctionType, class OperatorAssembler>
class BoxStokes
    : public NonlinearModel<Grid, Scalar, ProblemType, LocalJacobian, FunctionType, OperatorAssembler>
{
public:
    typedef Dune::NonlinearModel<Grid,
                                 Scalar,
                                 ProblemType,
                                 LocalJacobian,
                                 FunctionType,
                                 OperatorAssembler> NonlinearModel;

    BoxStokes(const Grid& grid, ProblemType& prob)
        : NonlinearModel(grid, prob), uOldTimeStep(grid)
    { }

    BoxStokes(const Grid& grid, ProblemType& prob, int level)
        : NonlinearModel(grid, prob, level), uOldTimeStep(grid, level)
    {     }

    virtual void initial() = 0;

    virtual void update(double& dt) = 0;

    virtual void solve() = 0;

    virtual ~BoxStokes () {}

    FunctionType uOldTimeStep;
};



/** \todo Please doc me! */

template<class Grid, class Scalar, int dim>
class LeafP1BoxStokes : public BoxStokes<Grid, Scalar, StokesProblem<Grid, Scalar>, BoxStokesJacobian<Grid, Scalar>,
                                         LeafP1Function<Grid, Scalar, dim+1>, LeafP1OperatorAssembler<Grid, Scalar, dim+1> >
{
public:
    enum{numEq = dim+1};

    typedef Grid GridType;

    // define the function type:
    typedef LeafP1Function<Grid, Scalar, numEq> FunctionType;

    // define the operator assembler type:
    typedef LeafP1OperatorAssembler<Grid, Scalar, numEq> OperatorAssembler;

    typedef Dune::BoxStokes<Grid, Scalar, StokesProblem<Grid, Scalar>, BoxStokesJacobian<Grid, Scalar>,
                            FunctionType, OperatorAssembler> BoxStokes;

    typedef LeafP1BoxStokes<Grid, Scalar, dim> ThisType;

    typedef BoxStokesJacobian<Grid, Scalar> LocalJacobian;

    // mapper: one data element per vertex
    template<int dimension>
    struct P1Layout
    {
        bool contains (Dune::GeometryType gt)
        {
            return gt.dim() == 0;
        }
    };

    typedef typename Grid::LeafGridView GV;
    typedef typename GV::IndexSet IS;
    typedef typename GV::IntersectionIterator IntersectionIterator;
    typedef MultipleCodimMultipleGeomTypeMapper<GV,P1Layout> VertexMapper;
    typedef typename ThisType::FunctionType::RepresentationType VectorType;
    typedef typename ThisType::OperatorAssembler::RepresentationType MatrixType;
    typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#ifdef HAVE_PARDISO
    //    SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
#endif

    LeafP1BoxStokes (const Grid& grid, StokesProblem<Grid, Scalar>& prob)
        : BoxStokes(grid, prob), grid_(grid), vertexmapper(grid.leafView()),
          size((*(this->u)).size()), pressure(size), xVelocity(size), yVelocity(size),
          uOldNewtonStep(size)
    { }

    VectorType& solOldNewtonStep()
    {
        return uOldNewtonStep;
    }

    void initial()
    {
        typedef typename Grid::Traits::template Codim<0>::Entity Element;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        enum{dimworld = Grid::dimensionworld};

        const GV& gridview(this->grid_.leafView());
        std::cout << "initializing solution." << std::endl;
        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
        {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Element& entity = *it;

            const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
                sfs=Dune::LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt, 1);
            int size = sfs.size();

            for (int i = 0; i < size; i++)
            {
                // get cell center in reference element
                const Dune::FieldVector<Scalar,dim>&local = sfs[i].position();

                // get global coordinate of cell center
                Dune::FieldVector<Scalar,dimworld> global = it->geometry().global(local);

                int globalId = vertexmapper.template map<dim>(entity, sfs[i].entity());

                for (int comp = 0; comp < dim; comp++)
                    (*(this->u))[globalId][comp] = 0;//this->problem.velocity(global)[comp];

                (*(this->u))[globalId][dim] = 0;
            }
        }
        //        (*(this->u))[116][dim] = 0.2;

        // set Dirichlet boundary conditions
        for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
        {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Element& entity = *it;

            const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type&
                sfs=Dune::LagrangeShapeFunctions<Scalar,Scalar,dim>::general(gt, 1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().fvGeom.update(entity);
            this->localJacobian().assembleBoundaryCondition(entity);

            IntersectionIterator endit = entity.ileafend();
            for (IntersectionIterator is = entity.ileafbegin();
                 is!=endit; ++is)
                if (is->boundary())
                {
                    for (int i = 0; i < size; i++)
                        // handle subentities of this face
                        for (int j = 0; j < ReferenceElements<Scalar,dim>::general(gt).size(is->indexInInside(), 1, sfs[i].codim()); j++)
                            if (sfs[i].entity() == ReferenceElements<Scalar,dim>::general(gt).subEntity(is->indexInInside(), 1, j, sfs[i].codim()))
                            {
                                if (this->localJacobian().bc(i)[1] == BoundaryConditions::dirichlet)
                                {
                                    // get cell center in reference element
                                    Dune::FieldVector<Scalar,dim> local = sfs[i].position();

                                    // get global coordinate of cell center
                                    Dune::FieldVector<Scalar,dimworld> global = it->geometry().global(local);

                                    int globalId = vertexmapper.template map<dim>(entity, sfs[i].entity());

                                    BoundaryConditions::Flags bctype = this->problem.bctype(global, entity, is, local);
                                    if (bctype == BoundaryConditions::dirichlet) {
                                        FieldVector<Scalar,numEq> dirichlet = this->problem.g(global, entity, is, local);
                                        //TODO: dim or numEq for the following loop??
                                        for (int eq = 0; eq < dim; eq++)
                                            (*(this->u))[globalId][eq] = dirichlet[eq];
                                    }
                                    else {
                                        std::cout << global << " is considered to be a Neumann node." << std::endl;
                                    }
                                }
                            }
                }
        }

        //        //set pressure condition (same node as in asseble())
        //        Iterator it = gridview.template begin<0>();
        //        unsigned int globalId = vertexmapper.template map<dim>(*it, 3);
        //
        //        //get global coordinates of node globalId
        //        // get geometry type
        //        Dune::GeometryType gt = it->geometry().type();
        //        // get local coordinates of node globalId (element it, node 3)
        //        Dune::FieldVector<Scalar,dim> local = ReferenceElements<Scalar, dim>::general(gt).position(3,dim);
        //        //std::cout<<"local:"<<local<<std::endl;
        //        // get global coordinates
        //        Dune::FieldVector<Scalar, dimworld> global = it->geometry().global(local);
        //        //std::cout<<"global:"<<global<<std::endl;
        //
        //        //set third position of u to pressure(globalID)
        //        //std::cout<<"pressure:"<<this->problem.pressure(global)<<std::endl;
        //        (*(this->u))[globalId][dim] = this->problem.pressure(global);
        //        //std::cout<<"pressure:"<<(*(this->u))[globalId][dim]<<std::endl;

        *(this->uOldTimeStep) = *(this->u);
        return;
    }


    void assemble()
    {
        *(this->f) = 0;
        this->localJacobian().clearVisited();
        this->A.assemble(this->localJacobian(), this->u, this->f);

        //        const GV& gridview(this->grid_.leafView());
        //        typedef typename GV::template Codim<0>::Iterator Iterator;
        //
        //        Iterator it = gridview.template begin<0>();
        //        unsigned int globalId = 81;//vertexmapper.template map<dim>(*it, 3);
        //
        //        MatrixType& A = *(this->A);
        //        for (typename MatrixType::RowIterator i=A.begin(); i!=A.end(); ++i)
        //            if(i.index()==globalId)
        //                for (typename MatrixType::ColIterator j=(*i).begin(); j!=(*i).end(); ++j)
        //                    A[i.index()][j.index()][dim] = 0.0;
        //        A[globalId][globalId][dim][dim] = 1.0;
        //        (*(this->f))[globalId][dim] = 0.0;
    }

    void update(double& dt)
    {
        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);
        double dtol = 1e-3;
        double rtol = 1e7;
        int maxIt = 10;
        double mindt = 1e-5;
        int goodIt = 4;
        int maxInc = 2;
        NewtonMethod<Grid, ThisType> newtonMethod(this->grid_, *this,
                                                  dtol, rtol, maxIt, mindt, goodIt, maxInc);
        newtonMethod.execute();
        dt = this->localJacobian().getDt();
        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    void solve()
    {
        // modify matrix and rhs for introducing pressure boundary condition
        // this is done here and not in the assembly step, since it is not needed for the coupled setting
        const GV& gridview(this->grid_.leafView());
        typedef typename GV::template Codim<0>::Iterator Iterator;

        Iterator it = gridview.template begin<0>();
        unsigned int globalId = vertexmapper.template map<dim>(*it, 3);

        MatrixType& A = *(this->A);
        for (typename MatrixType::RowIterator i=A.begin(); i!=A.end(); ++i)
            if(i.index()==globalId)
                for (typename MatrixType::ColIterator j=(*i).begin(); j!=(*i).end(); ++j)
                    A[i.index()][j.index()][dim] = 0.0;
        A[globalId][globalId][dim][dim] = 1.0;
        (*(this->f))[globalId][dim] = 0.0;

        Operator op(A);  // make operator out of matrix
        double red=1E-14;

#ifdef HAVE_PARDISO
        SeqPardiso<MatrixType,VectorType,VectorType> pardiso(A);
        LoopSolver<VectorType> solver(op, pardiso, red, 10, 2);
#else
        SeqILU0<MatrixType,VectorType,VectorType> ilu0(A,1.0);// a precondtioner
        BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,1);         // an inverse operator
#endif
        InverseOperatorResult r;
        solver.apply(*(this->u), *(this->f), r);

        return;
    }

    void globalDefect(FunctionType& defectGlobal) {
        typedef typename Grid::Traits::template Codim<0>::Entity Element;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename BoundaryConditions::Flags BCBlockType;

        const GV& gridview(this->grid_.leafView());
        (*defectGlobal)=0;

        // allocate flag vector to hold flags for essential boundary conditions
        std::vector<BCBlockType> essential(this->vertexmapper.size());
        for (typename std::vector<BCBlockType>::size_type i=0; i
                 <essential.size(); i++)
            essential[i] = BoundaryConditions::neumann;

        // iterate through leaf grid
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Element& entity = *it;
            this->localJacobian().fvGeom.update(entity);
            int size = this->localJacobian().fvGeom.numVertices;

            this->localJacobian().setLocalSolution(entity);
            this->localJacobian().computeElementData(entity);
            this->localJacobian().updateVariableData(entity, this->localJacobian().u);
            this->localJacobian().localDefect(entity, this->localJacobian().u);

            // begin loop over vertices
            for (int i=0; i < size; i++) {
                int globalId = this->vertexmapper.template map<dim>(entity,i);
                for (int equationnumber = 0; equationnumber < numEq; equationnumber++) {
                    if (this->localJacobian().bc(i)[equationnumber] == BoundaryConditions::neumann)
                        (*defectGlobal)[globalId][equationnumber]
                            += this->localJacobian().def[i][equationnumber];
                    else
                        essential[globalId] = BoundaryConditions::dirichlet;
                }
            }
        }

        for (typename std::vector<BCBlockType>::size_type i=0; i<essential.size(); i++)
            if (essential[i] == BoundaryConditions::dirichlet)
                for (int equationnumber = 0; equationnumber < dim; equationnumber++)
                    (*defectGlobal)[i][equationnumber] = 0;
    }

    void vtkout (const char* name, int k)
    {
        for (int i = 0; i < size; i++) {
            pressure[i] = (*(this->u))[i][dim];
            xVelocity[i] = (*(this->u))[i][0];
            yVelocity[i] = (*(this->u))[i][1];
        }

        VTKWriter<typename Grid::LeafGridView> vtkwriter(this->grid_.leafView());
        vtkwriter.addVertexData(pressure,"pressure");
        vtkwriter.addVertexData(xVelocity,"xVelocity");
        vtkwriter.addVertexData(yVelocity,"yVelocity");
        vtkwriter.write(name, VTKOptions::ascii);
    }

    const Grid& grid() const
    { return grid_; }



protected:
    const Grid& grid_;
    VertexMapper vertexmapper;
    int size;
    BlockVector<FieldVector<Scalar, 1> > pressure;
    BlockVector<FieldVector<Scalar, 1> > xVelocity;
    BlockVector<FieldVector<Scalar, 1> > yVelocity;
    VectorType uOldNewtonStep;
};

}
#endif
