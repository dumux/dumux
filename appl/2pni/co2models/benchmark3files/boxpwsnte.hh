// $Id$

#ifndef DUNE_BOXPWSNTE_HH
#define DUNE_BOXPWSNTE_HH

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
#include <dune/istl/paamg/amg.hh>
#include "dumux/pardiso/pardiso.hh"
#include "dumux/pardiso/identity.hh"
#include "dumux/2pni/2pnimodel.hh"
#include "dumux/2pni/2pniproblem.hh"
#include "dumux/2pni/fv/boxpwsntejacobian.hh"

#include "dumux/nonlinear/new_newtonmethod.hh"
#include "dumux/2pni/2pninewtoncontroller.hh"
#include "dumux/io/importfromdgf_leaf.hh"
#include <boost/format.hpp>

namespace Dune {

/**
   \brief Two phase model with Pw and Sn as primary unknowns

   This implements a two phase model with Pw and Sn as primary unknowns.
*/
template<class G, class RT, class VtkMultiWriter>
class BoxPwSnTe
    : public LeafP1TwoPhaseModel<G, RT, TwoPhaseHeatProblem<G, RT>, BoxPwSnTeJacobian<G, RT> >
{

public:
    // define the problem type (also change the template argument above)
    typedef TwoPhaseHeatProblem<G, RT> ProblemType;

    // define the local Jacobian (also change the template argument above)
    typedef BoxPwSnTeJacobian<G, RT> LocalJacobian;
    typedef LeafP1TwoPhaseModel<G, RT, ProblemType, LocalJacobian>
    ThisLeafP1TwoPhaseModel;
    typedef typename ThisLeafP1TwoPhaseModel::FunctionType FunctionType;

    typedef typename G::LeafGridView GV;
    typedef typename GV::IndexSet IS;

    enum {m = 3};

    typedef BoxPwSnTe<G, RT, VtkMultiWriter> ThisType;
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
        typedef RT                                                  Scalar;
        typedef typename ThisLeafP1TwoPhaseModel::FunctionType      Function;
        typedef typename ThisType::LocalJacobian                    LocalJacobian;
        typedef typename ThisLeafP1TwoPhaseModel::OperatorAssembler JacobianAssembler;
    };

    // HACK: traits for the domain of the problem. this is incomplete...
    struct DomainTraits {
        typedef RT   Scalar;
        typedef G    Grid;
    };

    typedef NewNewtonMethod<ThisType> NewtonMethod;
    typedef NewtonController<NewtonMethod> NewtonController;

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

    BoxPwSnTe(const G& g, ProblemType& prob)
        : ThisLeafP1TwoPhaseModel(g, prob)
    {}

    void initial() {
        typedef typename G::Traits::template Codim<0>::Entity Entity;
        typedef typename G::ctype DT;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename G::LeafGridView::IntersectionIterator IntersectionIterator;

        enum {dim = G::dimension};
        enum {dimworld = G::dimensionworld};

        const GV& gridview(this->_grid.leafView());

        this->localJacobian().outPressureN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outTemperature = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<RT, 1>(this->size);
        if(this->problem.soil().readPropertiesFlag() == true)
        {
            this->problem.soil().readSoilProperties();
            this->problem.soil().setSoilProperties();
        }
        this->localJacobian().outPermeabilityXDir = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outPorosity = vtkMultiWriter->template createField<RT, 1>(this->size);

        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Entity& entity = *it;

            this->localJacobian().fvGeom.update(entity);

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,1);
            int size = sfs.size();

            for (int i = 0; i < size; i++) {
                // get cell center in reference element
                const Dune::FieldVector<DT,dim>&local = sfs[i].position();

                // get global coordinate of cell center
                Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

                int globalId = this->vertexmapper.template map<dim>(entity,
                                                                    sfs[i].entity());

                // initialize cell concentration
                (*(this->u))[globalId] = this->problem.initial(global, entity, local);
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

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,
                                                                        1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().assembleBoundaryCondition(entity);

            IntersectionIterator endit = entity.ileafend();
            for (IntersectionIterator is = entity.ileafbegin(); is!=endit; ++is)
                if (is->boundary()) {
                    for (int i = 0; i < size; i++)
                        // handle subentities of this face
                        for (int j = 0; j < ReferenceElements<DT,dim>::general(gt).size(is->indexInInside(), 1, sfs[i].codim()); j++)
                            if (sfs[i].entity()
                                == ReferenceElements<DT,dim>::general(gt).subEntity(is->indexInInside(), 1,
                                                                                    j, sfs[i].codim())) {
                                for (int equationNumber = 0; equationNumber<m; equationNumber++) {
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
                                        FieldVector<int,m> dirichletIndex;
                                        FieldVector<BoundaryConditions::Flags, m>
                                            bctype = this->problem.bctype(
                                                                          global, entity, is,
                                                                          local);
                                        this->problem.dirichletIndex(global, entity, is,
                                                                     local, dirichletIndex);

                                        if (bctype[equationNumber]
                                            == BoundaryConditions::dirichlet) {
                                            FieldVector<RT,m>
                                                ghelp = this->problem.dirichlet(
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
        typedef typename G::Traits::template Codim<0>::Entity Entity;
        typedef typename G::ctype DT;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename G::LeafGridView::IntersectionIterator IntersectionIterator;

        enum {dim = G::dimension};
        enum {dimworld = G::dimensionworld};

        const GV& gridview(this->_grid.leafView());

        this->localJacobian().outPressureN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outTemperature = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<RT, 1>(this->size);
        if(this->problem.soil().readPropertiesFlag() == true)
        {
            this->problem.soil().readSoilProperties();
            this->problem.soil().setSoilProperties();
        }
        this->localJacobian().outPermeabilityXDir = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outPorosity = vtkMultiWriter->template createField<RT, 1>(this->size);

        int size = this->vertexmapper.size();
        Dune::BlockVector<FieldVector<double, m> > data(size);
        data=0;

        std::string restartFileName;
        restartFileName = (boost::format("data-%05d")
                           %restartNum).str();
        importFromDGF<GV>(data, restartFileName, false);

        for (int i=0;i<size;i++)
        {
            for (int j=0;j<m;j++)
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

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,1);
            int size = sfs.size();

            for (int i = 0; i < size; i++) {
                // get cell center in reference element
                const Dune::FieldVector<DT,dim>&local = sfs[i].position();

                // get global coordinate of cell center
                Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

                int globalId = this->vertexmapper.template map<dim>(entity,
                                                                    sfs[i].entity());
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

            const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<DT, RT, dim>::general(gt,
                                                                        1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().assembleBoundaryCondition(entity);

            IntersectionIterator endit = entity.ileafend();
            for (IntersectionIterator is = entity.ileafbegin(); is!=endit; ++is)
                if (is->boundary()) {
                    for (int i = 0; i < size; i++)
                        // handle subentities of this face
                        for (int j = 0; j < ReferenceElements<DT,dim>::general(gt).size(is->indexInInside(), 1, sfs[i].codim()); j++)
                            if (sfs[i].entity()
                                == ReferenceElements<DT,dim>::general(gt).subEntity(is->indexInInside(), 1,
                                                                                    j, sfs[i].codim())) {
                                for (int equationNumber = 0; equationNumber<m; equationNumber++) {
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
                                        FieldVector<int,m> dirichletIndex;
                                        FieldVector<BoundaryConditions::Flags, m>
                                            bctype = this->problem.bctype(
                                                                          global, entity, is,
                                                                          local);
                                        this->problem.dirichletIndex(global, entity, is,
                                                                     local, dirichletIndex);

                                        if (bctype[equationNumber]
                                            == BoundaryConditions::dirichlet) {
                                            FieldVector<RT,m>
                                                ghelp = this->problem.dirichlet(
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
        this->localJacobian().outPressureN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outTemperature = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outPermeabilityXDir = vtkMultiWriter->template createField<RT, 1>(this->size);
        this->localJacobian().outPorosity = vtkMultiWriter->template createField<RT, 1>(this->size);

        //        this->localJacobian().outPermeability = vtkMultiWriter->template createField<RT, 1>(this->size);

        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);

        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename G::Traits::template Codim<0>::Entity Entity;
        typedef typename G::ctype DT;
        enum {dim = G::dimension};

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
            }
            std::cout<<"timeStep resized to: "<<nextDt<<std::endl;
        }

        double mass(0);
        mass = this->totalCO2Mass();
        std::cout << mass<<"  /* total mass CO2 */ ";

        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    template<class MultiWriter>
    void addvtkfields (MultiWriter& writer)
    {
        //        BlockVector<FieldVector<RT, 1> > &xWN = *writer.template createField<RT, 1>(this->size);
        //        BlockVector<FieldVector<RT, 1> > &xAW = *writer.template createField<RT, 1>(this->size);
        //        BlockVector<FieldVector<RT, 1> > &satW = *writer.template createField<RT, 1>(this->size);

        //        writer.addScalarVertexFunction("nonwetting phase saturation", this->u, 1);
        writer.addScalarVertexFunction("pressure wetting phase", this->u, 0);
        writer.addVertexData(this->localJacobian().outPressureN,"pressure non-wetting phase");
        writer.addVertexData(this->localJacobian().outCapillaryP,"capillary pressure");
        writer.addVertexData(this->localJacobian().outTemperature,"temperature");
        writer.addVertexData(this->localJacobian().outSaturationW,"saturation wetting phase");
        writer.addVertexData(this->localJacobian().outSaturationN,"saturation non-wetting phase");
        writer.addVertexData(this->localJacobian().outDensityW,"density wetting phase");
        writer.addVertexData(this->localJacobian().outDensityN,"density non-wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityW,"mobility wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityN,"mobility non-wetting phase");
        writer.addVertexData(this->localJacobian().outPermeabilityXDir, "permeability in x direction");
        writer.addVertexData(this->localJacobian().outPorosity,"porosity");
    }

    void setVtkMultiWriter(VtkMultiWriter *writer)
    { vtkMultiWriter = writer; }

protected:
    VtkMultiWriter *vtkMultiWriter;
};

}
#endif
