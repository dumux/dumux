#ifndef DUNE_BOXSTOKES_HH
#define DUNE_BOXSTOKES_HH

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

namespace Dune
{
  template<class G, class RT, class ProblemType, class LocalJacobian,
            class FunctionType, class OperatorAssembler>
  class BoxStokes
  : public NonlinearModel<G, RT, ProblemType, LocalJacobian, FunctionType, OperatorAssembler>
  {
  public:
  typedef Dune::NonlinearModel<G,
                               RT,
                               ProblemType,
                               LocalJacobian,
                               FunctionType,
                               OperatorAssembler> NonlinearModel;

    BoxStokes(const G& g, ProblemType& prob)
    : NonlinearModel(g, prob), uOldTimeStep(g)
    { }

    BoxStokes(const G& g, ProblemType& prob, int level)
    : NonlinearModel(g, prob, level), uOldTimeStep(g, level)
    {     }

    virtual void initial() = 0;

    virtual void update(double& dt) = 0;

    virtual void solve() = 0;

    virtual ~BoxStokes () {}

    FunctionType uOldTimeStep;
  };





  template<class G, class RT, int dim>
  class LeafP1BoxStokes : public BoxStokes<G, RT, StokesProblem<G, RT>, BoxStokesJacobian<G, RT>,
                                        LeafP1FunctionExtended<G, RT, dim+1>, LeafP1OperatorAssembler<G, RT, dim+1> >
  {
  public:
      enum{numEq = dim+1};

      // define the function type:
      typedef LeafP1FunctionExtended<G, RT, numEq> FunctionType;

      // define the operator assembler type:
      typedef LeafP1OperatorAssembler<G, RT, numEq> OperatorAssembler;

      typedef Dune::BoxStokes<G, RT, StokesProblem<G, RT>, BoxStokesJacobian<G, RT>,
                              FunctionType, OperatorAssembler> BoxStokes;

      typedef LeafP1BoxStokes<G, RT, dim> ThisType;

      typedef BoxStokesJacobian<G, RT> LocalJacobian;

      // mapper: one data element per vertex
      template<int dimension>
      struct P1Layout
      {
          bool contains (Dune::GeometryType gt)
          {
              return gt.dim() == 0;
          }
      };

       typedef typename G::LeafGridView GV;
        typedef typename GV::IndexSet IS;
      typedef MultipleCodimMultipleGeomTypeMapper<G,IS,P1Layout> VertexMapper;
      typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator IntersectionIterator;
        typedef typename ThisType::FunctionType::RepresentationType VectorType;
        typedef typename ThisType::OperatorAssembler::RepresentationType MatrixType;
        typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#ifdef HAVE_PARDISO
//    SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
#endif

      LeafP1BoxStokes (const G& g, StokesProblem<G, RT>& prob)
      : BoxStokes(g, prob), grid_(g), vertexmapper(g, g.leafIndexSet()),
        size((*(this->u)).size()), pressure(size), xVelocity(size), yVelocity(size), count(0)
      { }

      virtual void initial()
      {
          typedef typename G::Traits::template Codim<0>::Entity Entity;
          typedef typename G::ctype DT;
          typedef typename GV::template Codim<0>::Iterator Iterator;
          enum{dimworld = G::dimensionworld};

          const GV& gridview(this->grid_.leafView());
          std::cout << "initializing solution." << std::endl;
          // iterate through leaf grid an evaluate c0 at cell center
          Iterator eendit = gridview.template end<0>();
          for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
          {
              // get geometry type
              Dune::GeometryType gt = it->geometry().type();

              // get entity
              const Entity& entity = *it;

              const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
                  sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt, 1);
              int size = sfs.size();

              IntersectionIterator is = IntersectionIteratorGetter<G,LeafTag>::end(entity);

              for (int i = 0; i < size; i++) {
                  int globalId = vertexmapper.template map<dim>(entity, sfs[i].entity());

              	// initialize cell concentration
                (*(this->u))[globalId] = 0;

/*                // get cell center in reference element
                Dune::FieldVector<DT,dim> local = sfs[i].position();

                // get global coordinate of cell center
                Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

                FieldVector<DT,dim> dirichlet = this->problem.g(global, entity, is, local);
            	for (int eq = 0; eq < dim; eq++)
            		(*(this->u))[globalId][eq] = dirichlet[eq];*/

              }
          }

          // set Dirichlet boundary conditions
          for (Iterator it = gridview.template begin<0>(); it != eendit; ++it)
          {
              // get geometry type
              Dune::GeometryType gt = it->geometry().type();

              // get entity
              const Entity& entity = *it;

              const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type&
                  sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt, 1);
              int size = sfs.size();

              // set type of boundary conditions
              this->localJacobian().fvGeom.update(entity);
              this->localJacobian().template assembleBC<LeafTag>(entity);

//               for (int i = 0; i < size; i++)
//                 std::cout << "bc[" << i << "] = " << this->localJacobian().bc(i) << std::endl;

              IntersectionIterator endit = IntersectionIteratorGetter<G,LeafTag>::end(entity);
              for (IntersectionIterator is = IntersectionIteratorGetter<G,LeafTag>::begin(entity);
                   is!=endit; ++is)
                  if (is->boundary())
                  {
                    for (int i = 0; i < size; i++)
                      // handle subentities of this face
                      for (int j = 0; j < ReferenceElements<DT,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
                        if (sfs[i].entity() == ReferenceElements<DT,dim>::general(gt).subEntity(is->numberInSelf(), 1, j, sfs[i].codim()))
                        {
                            if (this->localJacobian().bc(i)[1] == BoundaryConditions::dirichlet)
                            {
                                // get cell center in reference element
                                Dune::FieldVector<DT,dim> local = sfs[i].position();

                                // get global coordinate of cell center
                                Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

                                int globalId = vertexmapper.template map<dim>(entity, sfs[i].entity());

                                BoundaryConditions::Flags bctype = this->problem.bctype(global, entity, is, local);
//                                 std::cout << "global = " << global << ", id = " << globalId << std::endl;
                                if (bctype == BoundaryConditions::dirichlet) {
                                	FieldVector<DT,dim> dirichlet = this->problem.g(global, entity, is, local);
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

          *(this->uOldTimeStep) = *(this->u);
          //printvector(std::cout, *(this->u), "initial solution", "row", 200, 1, 3);
          return;
      }


    virtual void update(double& dt)
    {
        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);
        NewtonMethod<G, ThisType> newtonMethod(this->grid_, *this);
        newtonMethod.execute();
        dt = this->localJacobian().getDt();
        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    virtual void solve()
    {
    	count++;
    	MatrixType& A = *(this->A);
        //modify matrix for introducing pressure boundary condition
        const GV& gridview(this->grid_.leafView());
        typedef typename GV::template Codim<0>::Iterator Iterator;

        Iterator it = gridview.template begin<0>();
//        unsigned int globalId = 8;//vertexmapper.template map<dim>(*it, 3);
        unsigned int globalId = vertexmapper.template map<dim>(*it, 3);

        //std::cout << "globalId = " << globalId << std::endl;
    	for (typename MatrixType::RowIterator i=A.begin(); i!=A.end(); ++i)
            if(i.index()==globalId)
            	for (typename MatrixType::ColIterator j=(*i).begin(); j!=(*i).end(); ++j)
                       A[i.index()][j.index()][dim] = 0.0;
        (*(this->A))[globalId][globalId][dim][dim] = 1.0;
        (*(this->f))[globalId][dim] = 0.0;

        Operator op(A);  // make operator out of matrix
        double red=1E-18;

#ifdef HAVE_PARDISO
        //pardiso.factorize(A);
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
        typedef typename G::Traits::template Codim<0>::Entity Entity;
        typedef typename G::ctype DT;
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
            const Entity& entity = *it;
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
            	for (int equationnumber = 0; equationnumber < dim; equationnumber++)
            		(*defectGlobal)[i][equationnumber] = 0;
    }

    virtual void vtkout (const char* name, int k)
    {
        for (int i = 0; i < size; i++) {
            pressure[i] = (*(this->u))[i][dim];
            xVelocity[i] = (*(this->u))[i][0];
            yVelocity[i] = (*(this->u))[i][1];
        }

    	VTKWriter<typename G::LeafGridView> vtkwriter(this->grid_.leafView());
        vtkwriter.addVertexData(pressure,"pressure");
        vtkwriter.addVertexData(xVelocity,"xVelocity");
        vtkwriter.addVertexData(yVelocity,"yVelocity");
        vtkwriter.write(name, VTKOptions::ascii);
    }

    const G& grid() const
        { return grid_; }



protected:
  const G& grid_;
  VertexMapper vertexmapper;
  int size;
  BlockVector<FieldVector<RT, 1> > pressure;
  BlockVector<FieldVector<RT, 1> > xVelocity;
  BlockVector<FieldVector<RT, 1> > yVelocity;
  int count;
};

}
#endif
