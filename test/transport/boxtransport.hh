#ifndef DUNE_BOXTRANSPORT_HH
#define DUNE_BOXTRANSPORT_HH

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
#include "boxtransportjacobian.hh"
#include "dumux/transport/transportproblem.hh"
#include "dumux/pardiso/pardiso.hh"

namespace Dune
{
template<class G, class RT, class ProblemType, class LocalJacobian,
         class FunctionType, class OperatorAssembler>
class BoxTransport
    : public NonlinearModel<G, RT, ProblemType, LocalJacobian, FunctionType, OperatorAssembler>
{
public:
    typedef Dune::NonlinearModel<G,
                                 RT,
                                 ProblemType,
                                 LocalJacobian,
                                 FunctionType,
                                 OperatorAssembler> NonlinearModel;

    BoxTransport(const G& g, ProblemType& prob)
        : NonlinearModel(g, prob), uOldTimeStep(g)
    { }

    BoxTransport(const G& g, ProblemType& prob, int level)
        : NonlinearModel(g, prob, level), uOldTimeStep(g, level)
    {     }

    virtual void initial() = 0;

    virtual void update(double& dt) = 0;

    virtual void solve() = 0;

    FunctionType uOldTimeStep;
};





template<class G, class RT, class VC, int m=1>
class LeafP1BoxTransport : public BoxTransport<G, RT, TransportProblem<G, RT, VC>, BoxTransportJacobian<G, RT, VC>,
                                               LeafP1Function<G, RT, m>, LeafP1OperatorAssembler<G, RT, m> >
{
public:
    // define the function type:
    typedef LeafP1Function<G, RT, m> FunctionType;

    // define the operator assembler type:
    typedef LeafP1OperatorAssembler<G, RT, m> OperatorAssembler;

    typedef Dune::BoxTransport<G, RT, TransportProblem<G, RT, VC>, BoxTransportJacobian<G, RT, VC>,
                               FunctionType, OperatorAssembler> BoxTransport;

    typedef LeafP1BoxTransport<G, RT, VC, m> ThisType;

    typedef BoxTransportJacobian<G, RT, VC> LocalJacobian;

    // mapper: one data element per vertex
    template<int dim>
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
    SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
#endif

    LeafP1BoxTransport (const G& g, TransportProblem<G, RT, VC>& prob)
        : BoxTransport(g, prob), grid_(g), vertexmapper(g, g.leafIndexSet()),
          size((*(this->u)).size())
    { }

    virtual void initial()
    {
        typedef typename G::Traits::template Codim<0>::Entity Entity;
        typedef typename G::ctype DT;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        enum{dim = G::dimension};
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

                for (int i = 0; i < size; i++) {
                    // get cell center in reference element
                    const Dune::FieldVector<DT,dim>&local = sfs[i].position();

                    // get global coordinate of cell center
                    Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

                    int globalId = vertexmapper.template map<dim>(entity, sfs[i].entity());

                    // initialize cell concentration
                    (*(this->u))[globalId] = this->problem.initSat(global, entity, local);
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
                                            if (this->localJacobian().bc(i)[0] == BoundaryConditions::dirichlet)
                                                {
                                                    // get cell center in reference element
                                                    Dune::FieldVector<DT,dim> local = sfs[i].position();

                                                    // get global coordinate of cell center
                                                    Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

                                                    int globalId = vertexmapper.template map<dim>(entity, sfs[i].entity());

                                                    FieldVector<BoundaryConditions::Flags, m> bctype = this->problem.bctypeSat(global, entity, local);
                                                    //                                 std::cout << "global = " << global << ", id = " << globalId << std::endl;
                                                    if (bctype[0] == BoundaryConditions::dirichlet) {
                                                        (*(this->u))[globalId] = this->problem.dirichletSat(global, entity, local);
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
        NewtonMethod<G, ThisType> newtonMethod(this->grid_, *this);
        newtonMethod.execute();
        dt = this->localJacobian().getDt();
        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    virtual void solve()
    {
        Operator op(*(this->A));  // make operator out of matrix
        double red=1E-8;

#ifdef HAVE_PARDISO
        pardiso.factorize(*(this->A));
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
        typedef typename G::Traits::template Codim<0>::Entity Entity;
        typedef typename G::ctype DT;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        enum {dim = G::dimension};
        typedef array<BoundaryConditions::Flags, m> BCBlockType;

        const GV& gridview(this->grid_.leafView());
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
                for (int equationnumber = 0; equationnumber < m; equationnumber++) {
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
            for (int equationnumber = 0; equationnumber < m; equationnumber++) {
                if (essential[i][equationnumber] == BoundaryConditions::dirichlet)
                    (*defectGlobal)[i][equationnumber] = 0;
            }
    }

    virtual void vtkout (const char* name, int k)
    {
        VTKWriter<typename G::LeafGridView> vtkwriter(this->grid_.leafView());
        vtkwriter.addVertexData(*(this->u),"saturation");
        char fname[128];
        sprintf(fname, "%s-%05d", name, k);
        vtkwriter.write(fname, VTKOptions::ascii);
    }

    const G &grid() const
    { return grid_; }



protected:
    const G& grid_;
    VertexMapper vertexmapper;
    int size;
};

}
#endif
