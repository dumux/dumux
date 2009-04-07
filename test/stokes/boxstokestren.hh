#ifndef DUNE_BOXSTOKESTREN_HH
#define DUNE_BOXSTOKESTREN_HH

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
#include "boxstokesjacobian.hh"
#include "dumux/stokes/stokesproblem.hh"
#include "dumux/pardiso/pardiso.hh"

#include "boxstokestrenjacobian.hh"
#include "dumux/stokes/stokestrenproblem.hh"


namespace Dune
{
  template<class Grid, class Scalar, class ProblemType, class LocalJacobian,
            class FunctionType, class OperatorAssembler>
  class BoxStokesTrEn
  : public NonlinearModel<Grid, Scalar, ProblemType, LocalJacobian, FunctionType, OperatorAssembler>
  {
  public:
  typedef Dune::NonlinearModel<Grid,
                               Scalar,
                               ProblemType,
                               LocalJacobian,
                               FunctionType,
                               OperatorAssembler> NonlinearModel;

    BoxStokesTrEn(const Grid& grid, ProblemType& prob)
    : NonlinearModel(grid, prob), uOldTimeStep(grid)
    { }

    BoxStokesTrEn(const Grid& grid, ProblemType& prob, int level)
    : NonlinearModel(grid, prob, level), uOldTimeStep(grid, level)
    {     }

    virtual void initial() = 0;

    virtual void update(double& dt) = 0;

    virtual void solve() = 0;

    virtual ~BoxStokesTrEn () {}

    FunctionType uOldTimeStep;
  };

  template<class Grid, class Scalar, int dim>
  class LeafP1BoxStokesTrEn : public BoxStokesTrEn<Grid, Scalar, StokesTrEnProblem<Grid, Scalar>, BoxStokesTrEnJacobian<Grid, Scalar>,
                                        LeafP1Function<Grid, Scalar, dim+3>, LeafP1OperatorAssembler<Grid, Scalar, dim+3> >
  {
  public:
      enum{numEq = dim+3};

      typedef Grid GridType;

      // define the function type:
      typedef LeafP1Function<Grid, Scalar, numEq> FunctionType;

      // define the operator assembler type:
      typedef LeafP1OperatorAssembler<Grid, Scalar, numEq> OperatorAssembler;

      typedef Dune::BoxStokesTrEn<Grid, Scalar, StokesTrEnProblem<Grid, Scalar>, BoxStokesTrEnJacobian<Grid, Scalar>,
                              FunctionType, OperatorAssembler> BoxStokesTrEn;

      typedef LeafP1BoxStokesTrEn<Grid, Scalar, dim> ThisType;

      typedef BoxStokesTrEnJacobian<Grid, Scalar> LocalJacobian;

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
      typedef MultipleCodimMultipleGeomTypeMapper<Grid,IS,P1Layout> VertexMapper;
      typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;
      typedef typename ThisType::FunctionType::RepresentationType VectorType;
      typedef typename ThisType::OperatorAssembler::RepresentationType MatrixType;
      typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#ifdef HAVE_PARDISO
//    SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
#endif

      LeafP1BoxStokesTrEn (const Grid& grid, StokesTrEnProblem<Grid, Scalar>& prob)
      : BoxStokesTrEn(grid, prob), grid_(grid), vertexmapper(grid, grid.leafIndexSet()),
        size((*(this->u)).size()), pressure(size), xVelocity(size), yVelocity(size),
        partialDensity(size), temperature(size), uOldNewtonStep(size)
      { }

      VectorType& solOldNewtonStep()
      {
          return uOldNewtonStep;
      }

      virtual void initial()
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

              IntersectionIterator is = IntersectionIteratorGetter<Grid,LeafTag>::end(entity);

              for (int i = 0; i < size; i++)
              {
            	  // get cell center in reference element
            	  const Dune::FieldVector<Scalar,dim>&local = sfs[i].position();

            	  // get global coordinate of cell center
            	  Dune::FieldVector<Scalar,dimworld> global = it->geometry().global(local);

            	  int globalId = vertexmapper.template map<dim>(entity, sfs[i].entity());

            	  // initialize cell concentration
            	  (*(this->u))[globalId] = this->problem.initial(global, entity, local); //0;
              }
          }

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
              this->localJacobian().template assembleBC<LeafTag>(entity);

              IntersectionIterator endit = IntersectionIteratorGetter<Grid,LeafTag>::end(entity);
              for (IntersectionIterator is = IntersectionIteratorGetter<Grid,LeafTag>::begin(entity);
                   is!=endit; ++is)
                  if (is->boundary())
                  {
                    for (int i = 0; i < size; i++)
                      // handle subentities of this face
                      for (int j = 0; j < ReferenceElements<Scalar,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
                        if (sfs[i].entity() == ReferenceElements<Scalar,dim>::general(gt).subEntity(is->numberInSelf(), 1, j, sfs[i].codim()))
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
                                	FieldVector<Scalar,dim+2> dirichlet = this->problem.g(global, entity, is, local);
                                	for (int eq = 0; eq < dim+2; eq++)
                                		(*(this->u))[globalId][eq] = dirichlet[eq];
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


      virtual void assemble()
      {
          MatrixType& A = *(this->A);
          *(this->f) = 0;
          this->localJacobian().clearVisited();
          this->A.assemble(this->localJacobian(), this->u, this->f);

          //modify matrix for introducing pressure boundary condition
          const GV& gridview(this->grid_.leafView());
          typedef typename GV::template Codim<0>::Iterator Iterator;

          Iterator it = gridview.template begin<0>();
          unsigned int globalId = vertexmapper.template map<dim>(*it, 3);

          for (typename MatrixType::RowIterator i=A.begin(); i!=A.end(); ++i)
              if(i.index()==globalId)
                  for (typename MatrixType::ColIterator j=(*i).begin(); j!=(*i).end(); ++j)
                         A[i.index()][j.index()][dim+2] = 0.0;
          A[globalId][globalId][dim+2][dim+2] = 1.0;
          (*(this->f))[globalId][dim+2] = 0.0;
      }

    virtual void update(double& dt)
    {
        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);
        double dtol = 1e-10;
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

    virtual void solve()
    {
    	MatrixType& A = *(this->A);
        Operator op(A);  // make operator out of matrix
        double red=1E-18;

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

    virtual void globalDefect(FunctionType& defectGlobal) {
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
            this->localJacobian().template localDefect<LeafTag>(entity, this->localJacobian().u);

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
            	for (int equationnumber = 0; equationnumber < dim+2; equationnumber++)
            		(*defectGlobal)[i][equationnumber] = 0;
    }

    virtual void vtkout (const char* name, int k)
    {
        for (int i = 0; i < size; i++) {
            pressure[i] = (*(this->u))[i][dim+2];
            xVelocity[i] = (*(this->u))[i][0];
            yVelocity[i] = (*(this->u))[i][1];
            partialDensity[i]  = (*(this->u))[i][dim];
            temperature[i] = (*(this->u))[i][dim+1];
        }

    	VTKWriter<typename Grid::LeafGridView> vtkwriter(this->grid_.leafView());
        vtkwriter.addVertexData(pressure,"pressure");
        vtkwriter.addVertexData(xVelocity,"xVelocity");
        vtkwriter.addVertexData(yVelocity,"yVelocity");
        vtkwriter.addVertexData(partialDensity,"partialDensity");
        vtkwriter.addVertexData(temperature,"temperature");
        char fname[128];
        sprintf(fname, "%s-%05d", name, k);
        vtkwriter.write(fname, VTKOptions::ascii);
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
  BlockVector<FieldVector<Scalar, 1> > partialDensity;
  BlockVector<FieldVector<Scalar, 1> > temperature;
  VectorType uOldNewtonStep;
};

}
#endif
