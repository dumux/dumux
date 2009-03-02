// $Id$

#ifndef DUNE_BOXPNSW_HH
#define DUNE_BOXPNSW_HH

#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>

//#include <dune/common/array.hh>        // defines simple array class
#include <dune/common/fixedarray.hh>   // defines simple array classes
#include <dune/common/geometrytype.hh>
#include <dune/grid/sgrid.hh>          // a complete structured grid
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/common/universalmapper.hh>
#include <dune/grid/common/quadraturerules.hh>
#include <dune/common/collectivecommunication.hh>
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
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/disc/functions/functions.hh>
#include "dumux/operators/p1operatorextended.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/istl/paamg/amg.hh>
#include "dumux/pardiso/pardiso.hh"
#include "dumux/pardiso/identity.hh"
#include "dumux/twophase/twophasemodel.hh"
#include "dumux/twophase/twophaseproblem.hh"
#include "dumux/twophase/fv/boxpnswjacobian.hh"
#include "dumux/nonlinear/new_newtonmethod.hh"
#include "dumux/nonlinear/new_newtoncontroller.hh"
#include "dumux/io/importfromdgf_leaf.hh"
#include <boost/format.hpp>

namespace Dune {

/**
   \brief Two phase model with Pn and Sw as primary unknowns

   This implements a two phase model with Pn and Sw as primary unknowns.
*/
template<class GridT, class ScalarT, class VtkMultiWriter>
class BoxPnSw
    : public LeafP1TwoPhaseModel<GridT, ScalarT, TwoPhaseProblem<GridT, ScalarT>, BoxPnSwJacobian<GridT, ScalarT> >
{

public:

    typedef ScalarT Scalar;
    typedef GridT Grid;

    // define the problem type (also change the template argument above)
    typedef TwoPhaseProblem<Grid, Scalar> ProblemType;

    // define the local Jacobian (also change the template argument above)
    typedef BoxPnSwJacobian<Grid, Scalar> LocalJacobian;
    typedef Dune::LeafP1TwoPhaseModel<Grid, Scalar, ProblemType, LocalJacobian>
    ThisLeafP1TwoPhaseModel;
    typedef typename ThisLeafP1TwoPhaseModel::FunctionType FunctionType;

    typedef typename Grid::LeafGridView GV;
    typedef typename GV::IndexSet IS;

    enum {numEq = 2};

    typedef BoxPnSw<Grid, Scalar, VtkMultiWriter> ThisType;
    typedef typename ThisLeafP1TwoPhaseModel::FunctionType::RepresentationType VectorType;
    typedef typename ThisLeafP1TwoPhaseModel::OperatorAssembler::RepresentationType MatrixType;
    typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#if HAVE_MPI
#else
#ifdef HAVE_PARDISO
    SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
#endif
#endif

    //////////////////////
    // Stuff required for the new newton method

    //! The traits class for the new newton method.
    struct NewtonTraits {
        typedef ScalarT                                             Scalar;
        typedef typename ThisLeafP1TwoPhaseModel::FunctionType      Function;
        typedef typename ThisType::LocalJacobian                    LocalJacobian;
        typedef typename ThisLeafP1TwoPhaseModel::OperatorAssembler JacobianAssembler;
    };

    // HACK: traits for the domain of the problem. this is incomplete...
    struct DomainTraits {
        typedef ScalarT         Scalar;
        typedef GridT           Grid;
    };

    typedef Dune::NewNewtonMethod<ThisType> NewtonMethod;
    typedef Dune::NewtonController<NewtonMethod> NewtonController;

    typedef typename NewtonTraits::Function Function;
    Function &currentSolution()
    { return this->u; };

    LocalJacobian &getLocalJacobian()
    { return this->localJacobian(); }

    typedef typename NewtonTraits::JacobianAssembler JacobianAssembler;
    JacobianAssembler &jacobianAssembler()
    { return this->A; }
    // End of stuff for new newton method
    //////////////////////

    BoxPnSw(const Grid& grid, ProblemType& prob)
        : ThisLeafP1TwoPhaseModel(grid, prob)
    {}

    virtual void solve()
    {
        Operator op(*(this->A));  // make operator out of matrix
        double red=1E-14;

#if HAVE_MPI
        // set up parallel solvers
        typedef typename Grid::Traits::GlobalIdSet::IdType GlobalIdType;
        typedef OwnerOverlapCopyCommunication<GlobalIdType,int> CommunicationType;
        Dune::IndexInfoFromGrid<GlobalIdType,int> indexinfo;
        (this->u).fillIndexInfoFromGrid(indexinfo);
        typedef Dune::OwnerOverlapCopyCommunication<GlobalIdType,int> CommunicationType;
        CommunicationType oocc(indexinfo,MPIHelper::getCommunicator());
        int verbose=0;
        if (this->grid().comm().rank() == 0)
            verbose = 1;
        Dune::OverlappingSchwarzOperator<MatrixType,VectorType,VectorType,CommunicationType> oop(*(this->A),oocc);
        Dune::OverlappingSchwarzScalarProduct<VectorType,CommunicationType> osp(oocc);
        SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
        Dune::BlockPreconditioner<VectorType,VectorType,CommunicationType> parprec(ilu0,oocc);
        Dune::BiCGSTABSolver<VectorType> parcg(oop,osp,parprec,red,1000,verbose);

        // solve system
        Dune::InverseOperatorResult r;
        parcg.apply(*(this->u), *(this->f), r);
#else
#ifdef HAVE_PARDISO
        pardiso.factorize(*(this->A));
        BiCGSTABSolver<VectorType> solver(op,pardiso,red,100,2); // an inverse operator
        //    SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
        //LoopSolver<VectorType> solver(op, ilu0, red, 10, 2);
#else
        SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A), 1.0);// a precondtioner

        //SeqIdentity<MatrixType,VectorType,VectorType> ilu0(*(this->A));// a precondtioner
        BiCGSTABSolver<VectorType> solver(op, ilu0, red, 10000, 1); // an inverse operator
#endif
        InverseOperatorResult r;
        solver.apply(*(this->u), *(this->f), r);
#endif

        return;
    }

    virtual void initial()
    {
        typedef typename Grid::Traits::template Codim<0>::Entity Entity;
        typedef typename Grid::ctype DT;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

        enum {dim = Grid::dimension};
        enum {dimworld = Grid::dimensionworld};

        const GV& gridview(this->grid().leafView());

        // for multiwriter output
        this->localJacobian().outPressureW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);

        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Entity& entity = *it;

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, Scalar, dim>::general(gt,
                                                                            1);
            int size = sfs.size();

            for (int i = 0; i < size; i++) {
                // get cell center in reference element
                const Dune::FieldVector<DT,dim>&local = sfs[i].position();

                // get global coordinate of cell center
                Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

                int globalId = this->vertexmapper.template map<dim>(entity,
                                                                    sfs[i].entity());

                // initialize cell concentration
                (*(this->u))[globalId] = this->problem.initial(
                                                               global, entity, local);
            }
        }

        // set Dirichlet boundary conditions
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Entity& entity = *it;

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, Scalar, dim>::general(gt,
                                                                            1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().template assembleBC<LeafTag>(entity);

            IntersectionIterator
                endit = IntersectionIteratorGetter<Grid, LeafTag>::end(entity);
            for (IntersectionIterator is = IntersectionIteratorGetter<Grid,
                     LeafTag>::begin(entity); is!=endit; ++is)
                if (is->boundary()) {
                    for (int i = 0; i < size; i++)
                        // handle subentities of this face
                        for (int j = 0; j < ReferenceElements<DT,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
                            if (sfs[i].entity()
                                == ReferenceElements<DT,dim>::general(gt).subEntity(is->numberInSelf(), 1,
                                                                                    j, sfs[i].codim())) {
                                for (int equationNumber = 0; equationNumber<numEq; equationNumber++) {
                                    if (this->localJacobian().bc(i)[equationNumber]
                                        == BoundaryConditions::dirichlet) {
                                        // get cell center in reference element
                                        Dune::FieldVector<DT,dim>
                                            local = sfs[i].position();

                                        // get global coordinate of cell center
                                        Dune::FieldVector<DT,dimworld>
                                            global = it->geometry().global(local);

                                        int
                                            globalId = this->vertexmapper.template map<dim>(
                                                                                            entity, sfs[i].entity());
                                        FieldVector<int,numEq> dirichletIndex;
                                        FieldVector<BoundaryConditions::Flags, numEq>
                                            bctype = this->problem.bctype(
                                                                          global, entity, is,
                                                                          local);
                                        this->problem.dirichletIndex(global, entity, is,
                                                                     local, dirichletIndex);

                                        if (bctype[equationNumber]
                                            == BoundaryConditions::dirichlet) {
                                            FieldVector<Scalar,numEq>
                                                ghelp = this->problem.g(
                                                                        global, entity, is,
                                                                        local);
                                            (*(this->u))[globalId][dirichletIndex[equationNumber]]
                                                = ghelp[dirichletIndex[equationNumber]];
                                        }
                                    }
                                }
                            }
                }
            // to avoid error in paraview (otherwise initial time step cannot be displayed)
            this->localJacobian().setLocalSolution(entity);
            for (int i = 0; i < size; i++)
                this->localJacobian().updateVariableData(entity, this->localJacobian().u, i, false);
        }

        *(this->uOldTimeStep) = *(this->u);
        return;
    }

    virtual void restart(int restartNum=0)
    {
        typedef typename Grid::Traits::template Codim<0>::Entity Entity;
        typedef typename Grid::ctype DT;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

        enum {dim = Grid::dimension};
        enum {dimworld = Grid::dimensionworld};

        const GV& gridview(this->grid().leafView());

        // for multiwriter output
        this->localJacobian().outPressureW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);

        int size = this->vertexmapper.size();
        Dune::BlockVector<FieldVector<double, numEq> > data(size);
        data=0;

        // initialize primary variables
        std::string restartFileName;
        restartFileName = (boost::format("data-%05d")
                           %restartNum).str();
        importFromDGF<GV>(data, restartFileName, false);


        for (int i=0;i<size;i++)
        {
            for (int j=0;j<numEq;j++)
            {
                (*(this->u))[i][j]=data[i][j];
            }

        }

        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            //            const Entity& entity = *it;

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, Scalar, dim>::general(gt,
                                                                            1);
            int size = sfs.size();

            for (int i = 0; i < size; i++) {
                // get cell center in reference element
                const Dune::FieldVector<DT,dim>&local = sfs[i].position();

                // get global coordinate of cell center
                Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

                /*                int globalId = this->vertexmapper.template map<dim>(entity,
                                  sfs[i].entity());
                */

            }
        }

        // set Dirichlet boundary conditions
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Entity& entity = *it;

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,Scalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, Scalar, dim>::general(gt,
                                                                            1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().template assembleBC<LeafTag>(entity);

            IntersectionIterator
                endit = IntersectionIteratorGetter<Grid, LeafTag>::end(entity);
            for (IntersectionIterator is = IntersectionIteratorGetter<Grid,
                     LeafTag>::begin(entity); is!=endit; ++is)
                if (is->boundary()) {
                    for (int i = 0; i < size; i++)
                        // handle subentities of this face
                        for (int j = 0; j < ReferenceElements<DT,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
                            if (sfs[i].entity()
                                == ReferenceElements<DT,dim>::general(gt).subEntity(is->numberInSelf(), 1,
                                                                                    j, sfs[i].codim())) {
                                for (int equationNumber = 0; equationNumber<numEq; equationNumber++) {
                                    if (this->localJacobian().bc(i)[equationNumber]
                                        == BoundaryConditions::dirichlet) {
                                        // get cell center in reference element
                                        Dune::FieldVector<DT,dim>
                                            local = sfs[i].position();

                                        // get global coordinate of cell center
                                        Dune::FieldVector<DT,dimworld>
                                            global = it->geometry().global(local);

                                        int
                                            globalId = this->vertexmapper.template map<dim>(
                                                                                            entity, sfs[i].entity());
                                        FieldVector<int,numEq> dirichletIndex;
                                        FieldVector<BoundaryConditions::Flags, numEq>
                                            bctype = this->problem.bctype(
                                                                          global, entity, is,
                                                                          local);
                                        this->problem.dirichletIndex(global, entity, is,
                                                                     local, dirichletIndex);

                                        if (bctype[equationNumber]
                                            == BoundaryConditions::dirichlet) {
                                            FieldVector<Scalar,numEq>
                                                ghelp = this->problem.g(
                                                                        global, entity, is,
                                                                        local);
                                            (*(this->u))[globalId][dirichletIndex[equationNumber]]
                                                = ghelp[dirichletIndex[equationNumber]];
                                        }
                                    }
                                }
                            }
                }
            this->localJacobian().setLocalSolution(entity);
            for (int i = 0; i < size; i++)
                this->localJacobian().updateVariableData(entity, this->localJacobian().u, i, false);

        }

        *(this->uOldTimeStep) = *(this->u);
        return;
    }

    void updateModel (double& dt, double& nextDt)
    {
        this->localJacobian().outPressureW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        //        this->localJacobian.outPermeability = vtkMultiWriter->template createField<Scalar, 1>(this->size);

        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);

        // execute newton method
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename Grid::Traits::template Codim<0>::Entity Entity;
        typedef typename Grid::ctype DT;
        enum {dim = Grid::dimension};

        bool newtonLoop = false;
        while(!newtonLoop)
        {
            nextDt = this->localJacobian().getDt();
            NewtonMethod newton(*this); // *this means object itself (boxpwsn)
            NewtonController newtonCtl(1e-7, 6);
            newtonLoop = newton.execute(*this, newtonCtl);
            nextDt = newtonCtl.suggestTimeStepSize(nextDt);
            this->localJacobian().setDt(nextDt);
            if(!newtonLoop){
                *this->u = *this->uOldTimeStep;
            }
            std::cout<<"timeStep resized to: "<<nextDt<<std::endl;
        }

        //        double Flux(0), Mass(0);
        //        Flux = this->computeFlux();
        //        Mass = totalCO2Mass();

        this->localJacobian().clearVisited();
        //        std::cout << Flux << ", "<< Mass;

        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    template<class MultiWriter>
    void addvtkfields (MultiWriter& writer)
    {
        writer.addScalarVertexFunction("pressure non-wetting phase", this->u, 0);
        writer.addVertexData(this->localJacobian().outPressureW,"pressure wetting phase");
        writer.addVertexData(this->localJacobian().outCapillaryP,"capillary pressure");
        writer.addVertexData(this->localJacobian().outSaturationW,"saturation wetting phase");
        writer.addVertexData(this->localJacobian().outSaturationN,"saturation non-wetting phase");
        writer.addVertexData(this->localJacobian().outDensityW,"density wetting phase");
        writer.addVertexData(this->localJacobian().outDensityN,"density non-wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityW,"mobility wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityN,"mobility non-wetting phase");
    }

    void setVtkMultiWriter(VtkMultiWriter *writer)
    { vtkMultiWriter = writer; }

protected:
    VtkMultiWriter *vtkMultiWriter;
};

}
#endif


