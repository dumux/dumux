#ifndef DUNE_PARALLELBOXDIFFUSION_HH
#define DUNE_PARALLELBOXDIFFUSION_HH

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
#include <dune/istl/schwarz.hh>
#include <dumux/operators/owneroverlapcopyextended.hh>
#include "dumux/nonlinear/nonlinearmodel.hh"
#include "dumux/pardiso/pardiso.hh"
#include "dumux/fvgeometry/fvelementgeometry.hh"
#include "dumux/nonlinear/newtonmethod.hh"
#include "parallelboxdiffusionjacobian.hh"
#include "diffusionparameters.hh"
#include "schwarzcommunication.hh"
#include "schwarzpreconditioner.hh"
#include "dumux/pardiso/pardiso.hh"
#include "dumux/operators/p1operatorextended.hh"

namespace Dune
{
template<class G, class RT, class ProblemType, class LocalJacobian,
class FunctionType, class OperatorAssembler>
class ParallelBoxDiffusion
: public NonlinearModel<G, RT, ProblemType, LocalJacobian, FunctionType, OperatorAssembler>
{
public:
    typedef NonlinearModel<G, RT, ProblemType, LocalJacobian,
    FunctionType, OperatorAssembler> NonlinearModel;

    ParallelBoxDiffusion(const G& g, ProblemType& prob)
    : NonlinearModel(g, prob), uOldTimeStep(g)
    { }

    ParallelBoxDiffusion(const G& g, ProblemType& prob, int level)
    : NonlinearModel(g, prob, level), uOldTimeStep(g, level)
    {     }

    virtual void initial() = 0;

    virtual void update(double& dt) = 0;

    virtual void solve() = 0;

    FunctionType uOldTimeStep;
};





template<class G, class RT, int m=1>
class LeafP1ParallelBoxDiffusion : public ParallelBoxDiffusion<G, RT, DiffusionParameters<G, RT>, ParallelBoxDiffusionJacobian<G, RT>,
LeafP1FunctionExtended<G, RT, m>, LeafP1OperatorAssembler<G, RT, m> >
{
public:
    // define the function type:
        typedef LeafP1FunctionExtended<G, RT, m> FunctionType;

        // define the operator assembler type:
        typedef LeafP1OperatorAssembler<G, RT, m> OperatorAssembler;

        typedef ParallelBoxDiffusion<G, RT, DiffusionParameters<G, RT>, ParallelBoxDiffusionJacobian<G, RT>,
        FunctionType, OperatorAssembler> ParallelBoxDiffusion;

        typedef LeafP1ParallelBoxDiffusion<G, RT, m> ThisType;

        typedef ParallelBoxDiffusionJacobian<G, RT> LocalJacobian;

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
        typedef typename G::Traits::GlobalIdSet::IdType GlobalIdType;
        typedef OwnerOverlapCopyExtendedCommunication<GlobalIdType,int> CommunicationType;
//#ifdef HAVE_PARDISO
//    SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
//#endif

        LeafP1ParallelBoxDiffusion (const G& g, DiffusionParameters<G, RT>& prob)
        : ParallelBoxDiffusion(g, prob), grid(g), vertexmapper(g, g.leafIndexSet())
        { }

        virtual void initial()
        {
            typedef typename G::Traits::template Codim<0>::Entity Entity;
            typedef typename G::ctype DT;
            typedef typename GV::template Codim<0>::Iterator Iterator;
            enum{dim = G::dimension};
            enum{dimworld = G::dimensionworld};

            *(this->u) = 0;

            const GV& gridview(this->grid.leafView());

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
                    int globalId = vertexmapper.template map<dim>(entity, sfs[i].entity());

                    // initialize cell concentration
                    (*(this->u))[globalId] = 0;
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
                this->localJacobian.fvGeom.update(entity);
                this->localJacobian.template assembleBC<LeafTag>(entity);

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
                                    if (this->localJacobian.bc(i)[0] == BoundaryConditions::dirichlet)
                                    {
                                        // get cell center in reference element
                                        Dune::FieldVector<DT,dim> local = sfs[i].position();

                                        // get global coordinate of cell center
                                        Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

                                        int globalId = vertexmapper.template map<dim>(entity, sfs[i].entity());

                                        FieldVector<BoundaryConditions::Flags, m> bctype = this->problem.bctype(global, entity, is, local);

                                        if (bctype[0] == BoundaryConditions::dirichlet)
                                            (*(this->u))[globalId] = this->problem.g(global, entity, is, local);
                                    }
                                }
                    }
            }

            *(this->uOldTimeStep) = *(this->u);
            return;
        }


        virtual void update(double& dt)
        {
            this->localJacobian.setDt(dt);
            this->localJacobian.setOldSolution(this->uOldTimeStep);
            NewtonMethod<G, ThisType> newtonMethod(this->grid, *this);
            newtonMethod.execute();
            dt = this->localJacobian.getDt();
            *(this->uOldTimeStep) = *(this->u);

            return;
        }

        virtual void solve()
        {
            typedef typename G::Traits::GlobalIdSet::IdType GlobalIdType;
            typedef typename Dune::LeafP1FunctionExtended<G,RT>::P1IndexInfoFromGrid P1IndexInfoFromGrid;

            Dune::MatrixAdapter<MatrixType,VectorType,VectorType> op(*(this->A));
            //SeqPardiso<MatrixType,VectorType,VectorType> ilu0;
            //ilu0.factorize(*(this->A));
//#ifdef HAVE_PARDISO
//            pardiso.factorize(*(this->A));
//#else
            Dune::SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);
//#endif
            //Dune::SeqSSOR<MatrixType,VectorType,VectorType> ssor(*(this->A),1,0.8);

#if HAVE_MPI
            // set up parallel solvers
            Dune::IndexInfoFromGrid<GlobalIdType,int> indexinfo;
            (this->u).fillIndexInfoFromGrid(indexinfo);
            typedef Dune::OwnerOverlapCopyExtendedCommunication<GlobalIdType,int> CommunicationType;
            CommunicationType oocc(indexinfo,grid.comm());
            int verbose=0;
            if (grid.comm().rank() == 0)
                verbose = 1;
            Dune::OverlappingSchwarzOperator<MatrixType,VectorType,VectorType,CommunicationType> oop(*(this->A),oocc);
            Dune::OverlappingSchwarzScalarProduct<VectorType,CommunicationType> osp(oocc);
//#ifdef HAVE_PARDISO
//            Dune::BlockPreconditioner<VectorType,VectorType,CommunicationType> parprec(pardiso,oocc);
//#else
            Dune::BlockPreconditioner<VectorType,VectorType,CommunicationType> parprec(ilu0,oocc);
//#endif
            //Dune::LoopSolver<VectorType> parcg(oop,osp,parprec,1E-12,1000,verbose);
            Dune::CGSolver<VectorType> parcg(oop,osp,parprec,1E-20, 1000,verbose);

            // solve system
            Dune::InverseOperatorResult r;
            parcg.apply(*(this->u), *(this->f), r);
#endif
            return;
        }

        virtual void vtkout (const char* name, int k)
        {
//            VTKWriter<typename G::LeafGridView> vtkwriter(this->grid.leafView());
//            vtkwriter.addVertexData(*(this->u),"pressure");
//            vtkwriter.write(name, VTKOptions::ascii);
        }

        virtual void globalDefect(FunctionType& defectGlobal) {
            typedef typename G::Traits::template Codim<0>::Entity Entity;
            typedef typename G::ctype DT;
            typedef typename GV::template Codim<0>::Iterator Iterator;
            enum {dim = G::dimension};
            typedef array<BoundaryConditions::Flags, m> BCBlockType;

            const GV& gridview(this->grid.leafView());
            (*defectGlobal)=0;

#if HAVE_MPI
            IndexInfoFromGrid<GlobalIdType,int> indexinfo;
            (this->u).fillIndexInfoFromGrid(indexinfo);
            CommunicationType oocc(indexinfo,grid.comm());
#endif

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
                this->localJacobian.fvGeom.update(entity);
                int size = this->localJacobian.fvGeom.numVertices;

                this->localJacobian.setLocalSolution(entity);
                this->localJacobian.computeElementData(entity);
                this->localJacobian.updateVariableData(entity, this->localJacobian.u);
                this->localJacobian.template localDefect<LeafTag>(entity, this->localJacobian.u);

                // begin loop over vertices
                for (int i=0; i < size; i++) {
                    int globalId = this->vertexmapper.template map<dim>(entity,i);
                    for (int equationnumber = 0; equationnumber < m; equationnumber++) {
                        if (this->localJacobian.bc(i)[equationnumber] == BoundaryConditions::neumann)
                            (*defectGlobal)[globalId][equationnumber]
                                    += this->localJacobian.def[i][equationnumber];
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

            //oocc.addAllToAll(*defectGlobal, *defectGlobal);

            //std::cout << grid.comm().rank() << ": norm(defect) = " << oocc.norm(*defectGlobal) << std::endl;

        }

        virtual double residual(FunctionType& defectGlobal)
        {
            globalDefect(defectGlobal);

#if HAVE_MPI
            IndexInfoFromGrid<GlobalIdType,int> indexinfo;
            (this->u).fillIndexInfoFromGrid(indexinfo);
            CommunicationType oocc(indexinfo,grid.comm());
            return oocc.norm(*defectGlobal);
#endif

            return 1e0;
        }

protected:
    const G& grid;
    VertexMapper vertexmapper;
};

}
#endif
