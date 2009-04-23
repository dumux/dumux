// $Id: boxco2.hh 766 2008-11-05 13:10:26Z melanie $

#ifndef DUNE_BOXCO2_HH
#define DUNE_BOXCO2_HH

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
#include <dune/grid/common/scsgmapper.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/disc/functions/functions.hh>
#include <dune/disc/operators/p1operator.hh>
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/paamg/amg.hh>
#include "dumux/pardiso/pardiso.hh"
#include "dumux/minc2p2cni/minc2p2cnimodel.hh"
#include "dumux/minc2p2cni/minc2p2cniproblem.hh"
#include "dumux/minc2p2cni/fv/boxmincco2jacobian.hh"

#include "dumux/nonlinear/new_newtonmethod.hh"
#include "dumux/minc2p2cni/minc2p2cninewtoncontroller.hh"
#include "dumux/io/importfromdgf_leaf.hh"
namespace Dune
{
/**
   \brief Non-isothermal two phase two component model with Pw, Sn/X and Temp as primary unknowns

   This implements a non-isothermal two phase two component model with Pw, Sn/X and Temp as primary unknowns
*/
template<class ThisGrid, class ThisScalar, int numEq, class VtkMultiWriter>
class BoxMincCO2
    : public LeafP1TwoPhaseModel<ThisGrid, ThisScalar, TwoPTwoCNIProblem<ThisGrid, ThisScalar, numEq>,
                                 BoxMincCO2Jacobian<ThisGrid, ThisScalar, numEq>, numEq >
{
public:
    // define the problem type (also change the template argument above)
    typedef TwoPTwoCNIProblem<ThisGrid, ThisScalar, numEq> ProblemType;

    // define the local Jacobian (also change the template argument above)
    typedef BoxMincCO2Jacobian<ThisGrid, ThisScalar, numEq> LocalJacobian;
    typedef LeafP1TwoPhaseModel<ThisGrid, ThisScalar, ProblemType, LocalJacobian, numEq> ThisLeafP1TwoPhaseModel;
    typedef BoxMincCO2<ThisGrid, ThisScalar, numEq, VtkMultiWriter> ThisType;

    typedef typename ThisLeafP1TwoPhaseModel::FunctionType FunctionType;

    typedef typename ThisGrid::LeafGridView GridView;

    typedef typename ThisLeafP1TwoPhaseModel::FunctionType::RepresentationType VectorType;
    typedef typename ThisLeafP1TwoPhaseModel::OperatorAssembler::RepresentationType MatrixType;
    typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#ifdef HAVE_PARDISO
    SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
#endif

    //////////////////////
    // Stuff required for the new newton method

    //! The traits class for the new newton method.
    struct NewtonTraits {
        typedef ThisScalar                                                  Scalar;
        typedef typename ThisLeafP1TwoPhaseModel::FunctionType      Function;
        typedef typename ThisType::LocalJacobian                    LocalJacobian;
        typedef typename ThisLeafP1TwoPhaseModel::OperatorAssembler JacobianAssembler;
    };

    // HACK: traits for the domain of the problem. this is incomplete...
    struct DomainTraits {
        typedef ThisScalar   Scalar;
        typedef ThisGrid Grid;
    };

    typedef NewNewtonMethod<ThisType> NewtonMethod;
    typedef TwoPTwoCNINewtonController<NewtonMethod> NewtonController;

    typedef typename NewtonTraits::Function Function;
    Function &currentSolution()
    { return this->u; };

    LocalJacobian &getLocalJacobian()
    { return this->localJacobian; }

    typedef typename NewtonTraits::JacobianAssembler JacobianAssembler;
    JacobianAssembler &jacobianAssembler()
    { return this->A; }
    // End of stuff for new newton method
    //////////////////////


    BoxMincCO2(const ThisGrid& g, ProblemType& prob)
        : ThisLeafP1TwoPhaseModel(g, prob)// (this->size) vectors
    { }

    void initial() {
        typedef typename ThisGrid::Traits::template Codim<0>::Entity Entity;
        typedef typename ThisGrid::ctype DT;
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        typedef typename ThisGrid::LeafGridView::IntersectionIterator IntersectionIterator;

        enum {dim = ThisGrid::dimension};
        enum {dimworld = ThisGrid::dimensionworld};

        const GridView& gridview(this->_grid.leafView());

        this->localJacobian().outPressureN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outTemperature = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracAir = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracWater = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);

        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Entity& entity = *it;

            this->localJacobian().fvGeom.update(entity);

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,ThisScalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, ThisScalar, dim>::general(gt,1);
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

                // initialize variable phaseState
                this->localJacobian().sNDat[globalId].phaseState =
                    this->problem.initialPhaseState(global, entity, local);
                // initialize variable oldPhaseState
                this->localJacobian().sNDat[globalId].oldPhaseState =
                    this->problem.initialPhaseState(global, entity, local);

            }
            this->localJacobian().clearVisited();
            this->localJacobian().initiateStaticData(entity);
        }

        // set Dirichlet boundary conditions
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Entity& entity = *it;

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,ThisScalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, ThisScalar, dim>::general(gt,
                                                                                1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().assembleBoundaryCondition(entity);

            IntersectionIterator
                endit = IntersectionIteratorGetter<ThisGrid, LeafTag>::end(entity);
            for (IntersectionIterator is = IntersectionIteratorGetter<ThisGrid,
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
                                            FieldVector<ThisScalar,numEq>
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

    void restart(int restartNum=0) {
        typedef typename ThisGrid::Traits::template Codim<0>::Entity Entity;
        typedef typename ThisGrid::ctype DT;
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        typedef typename ThisGrid::LeafGridView::IntersectionIterator IntersectionIterator;

        enum {dim = ThisGrid::dimension};
        enum {dimworld = ThisGrid::dimensionworld};

        const GridView& gridview(this->_grid.leafView());

        this->localJacobian().outPressureN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outTemperature = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracAir = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracWater = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);

        int size = this->vertexmapper.size();
        Dune::BlockVector<FieldVector<double, numEq+1> > data(size);
        data=0;

        importFromDGF<GridView>(data, "data", false);

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
            const Entity& entity = *it;

            this->localJacobian().fvGeom.update(entity);

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,ThisScalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, ThisScalar, dim>::general(gt,1);
            int size = sfs.size();

            for (int i = 0; i < size; i++) {
                // get cell center in reference element
                const Dune::FieldVector<DT,dim>&local = sfs[i].position();

                // get global coordinate of cell center
                Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

                int globalId = this->vertexmapper.template map<dim>(entity,
                                                                    sfs[i].entity());

                // initialize cell concentration


                // initialize variable phaseState
                this->localJacobian().sNDat[globalId].phaseStateF =
                    data[globalId][numEq];
                // initialize variable oldPhaseState
                this->localJacobian().sNDat[globalId].oldPhaseStateF = data[globalId][numEq];

            }
            this->localJacobian().clearVisited();
            this->localJacobian().initiateStaticData(entity);
        }

        // set Dirichlet boundary conditions
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Entity& entity = *it;

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,ThisScalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, ThisScalar, dim>::general(gt,
                                                                                1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().assembleBoundaryCondition(entity);

            IntersectionIterator
                endit = IntersectionIteratorGetter<ThisGrid, LeafTag>::end(entity);
            for (IntersectionIterator is = IntersectionIteratorGetter<ThisGrid,
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
                                            FieldVector<ThisScalar,numEq>
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
    void updateState()
    {
        typedef typename ThisGrid::Traits::template Codim<0>::Entity Entity;
        typedef typename ThisGrid::ctype DT;
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        typedef typename ThisGrid::LeafGridView::IntersectionIterator IntersectionIterator;

        enum {dim = ThisGrid::dimension};
        enum {dimworld = ThisGrid::dimensionworld};

        const GridView& gridview(this->_grid.leafView());
        // iterate through leaf grid and evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {

            const Entity& entity = *it;
            this->localJacobian().fvGeom.update(entity);
            this->localJacobian().setLocalSolution(entity);
            this->localJacobian().computeElementData(entity);
            this->localJacobian().updateStaticData(entity, this->localJacobian().u);
        }
        return;
    }

    virtual void globalDefect(FunctionType& defectGlobal)
    {
        ThisLeafP1TwoPhaseModel::globalDefect(defectGlobal);
    }

    void solve()
    {
        Operator op(*(this->A));  // make operator out of matrix
        double red=1E-9;

#ifdef HAVE_PARDISO
        //    SeqPardiso<MatrixType,VectorType,VectorType> ilu0(*(this->A));
        pardiso.factorize(*(this->A));
        BiCGSTABSolver<VectorType> solver(op,pardiso,red,100,2);         // an inverse operator
        //    SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
        //LoopSolver<VectorType> solver(op, ilu0, red, 10, 2);
#else
        SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner

        //SeqIdentity<MatrixType,VectorType,VectorType> ilu0(*(this->A));// a precondtioner
        BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,1);         // an inverse operator
#endif
        InverseOperatorResult r;
        solver.apply(*(this->u), *(this->f), r);

        return;
    }


    void update(double& dt)
    {
        DUNE_THROW(NotImplemented, "the update method is deprecated. use updateModel()");
    }


    void updateModel(double& dt, double &nextDt)
    {
        this->localJacobian().outPressureN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outTemperature = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracAir = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMassFracWater = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);

        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);

        typedef typename GridView::template Codim<0>::Iterator Iterator;
        typedef typename ThisGrid::Traits::template Codim<0>::Entity Entity;
        typedef typename ThisGrid::ctype DT;
        enum {dim = ThisGrid::dimension};

        ///////////////////////////////////
        // define solver tolerances here
        ///////////////////////////////////

        /*
//////////////
ThisScalar absTol = 1e-1;
ThisScalar relTol = 1e-9;
///////////////////
NewtonMethod<ThisGrid, ThisType> newtonMethod(this->_grid, *this, relTol, absTol);
newtonMethod.execute();
        */
        bool newtonLoop = false;
        while(!newtonLoop)
        {
            nextDt = this->localJacobian().getDt();
            NewtonMethod newton(*this);
            NewtonController newtonCtl;
            newtonLoop = newton.execute(*this, newtonCtl);
            nextDt = newtonCtl.suggestTimeStepSize(nextDt);
            this->localJacobian().setDt(nextDt);
            if(!newtonLoop){
                *this->u = *this->uOldTimeStep;
                this->localJacobian().resetPhaseState();
            }
            std::cout<<"timeStep resized to: "<<nextDt<<std::endl;
        }

        double Flux(0), Mass(0);
        Flux = this->computeFlux();
        Mass = this->totalCO2Mass();

        this->localJacobian().updatePhaseState(); // update variable oldPhaseState
        this->localJacobian().clearVisited();
        updateState();                            // phase switch after each timestep
        std::cout << Flux << ", "<< Mass;

        *(this->uOldTimeStep) = *(this->u);

        // update old phase state for computation of ComputeM(..uold..)

        return;
    }


    template<class MultiWriter>
    void addvtkfields (MultiWriter& writer)
    {
        //        BlockVector<FieldVector<ThisScalar, 1> > &xWN = *writer.template createField<ThisScalar, 1>(this->size);
        //        BlockVector<FieldVector<ThisScalar, 1> > &xAW = *writer.template createField<ThisScalar, 1>(this->size);
        //        BlockVector<FieldVector<ThisScalar, 1> > &satW = *writer.template createField<ThisScalar, 1>(this->size);

        //        writer.addScalarVertexFunction("nonwetting phase saturation", this->u, 1);
        writer.addScalarVertexFunction("pressure wetting phase", this->u, 0);

        ////        writer.addScalarVertexFunction("nonwetting phase saturation",
        ////                                        this->u,
        ////                                        1);
        //        writer.addScalarVertexFunction("wetting phase pressure",
        //                                        this->u,
        //                                        0);
        //        writer.addVertexData(&satW,"wetting phase saturation");
        writer.addVertexData(this->localJacobian().outPressureN,"pressure non-wetting phase");
        writer.addVertexData(this->localJacobian().outCapillaryP,"capillary pressure");
        writer.addVertexData(this->localJacobian().outTemperature,"temperature");
        writer.addVertexData(this->localJacobian().outSaturationW,"saturation wetting phase");
        writer.addVertexData(this->localJacobian().outSaturationN,"saturation non-wetting phase");
        writer.addVertexData(this->localJacobian().outMassFracAir,"massfraction co2 in wetting phase");
        writer.addVertexData(this->localJacobian().outMassFracWater,"massfraction water in non-wetting phase");
        writer.addVertexData(this->localJacobian().outDensityW,"density wetting phase");
        writer.addVertexData(this->localJacobian().outDensityN,"density non-wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityW,"mobility wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityN,"mobility non-wetting phase");
        writer.addVertexData(this->localJacobian().outPhaseState,"phase state");

        //        writer.addVertexData(&xWN, "water in co2");
        //        writer.addVertexData(&xAW, "dissolved co2");
    }



    void vtkout (const char* name, int k)
    {

    }


    void setVtkMultiWriter(VtkMultiWriter *writer)
    { vtkMultiWriter = writer; }
protected:
    VtkMultiWriter *vtkMultiWriter;

};
}
#endif
