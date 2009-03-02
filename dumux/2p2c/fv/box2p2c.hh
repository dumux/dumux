// $Id$

#ifndef DUNE_BOX2P2C_HH
#define DUNE_BOX2P2C_HH

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
#include "dumux/operators/p1operatorextended.hh"
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/paamg/amg.hh>
#include "dumux/pardiso/pardiso.hh"
#include "dumux/twophase/twophasemodel.hh"
#include "dumux/2p2c/2p2cproblem.hh"
#include "dumux/2p2c/fv/box2p2cjacobian.hh"

#include "dumux/nonlinear/new_newtonmethod.hh"
#include "dumux/2p2c/2p2cnewtoncontroller.hh"
#include "dumux/io/importfromdgf_leaf.hh"

namespace Dune
{
/**
   \brief Isothermal two phase two component model with Pw and Sn/X as primary unknowns

   This implements an isothermal two phase two component model with Pw and Sn/X as primary unknowns
*/
template<class GridT, class ScalarT, class VtkMultiWriter>
class Box2P2C
    : public LeafP1TwoPhaseModel<GridT, ScalarT, TwoPTwoCProblem<GridT, ScalarT>, Box2P2CJacobian<GridT, ScalarT> >
{
public:
    // define the problem type (also change the template argument above)
    typedef ScalarT Scalar;
    typedef GridT   Grid;
    typedef TwoPTwoCProblem<Grid, Scalar> ProblemType;

    // define the local Jacobian (also change the template argument above)
    typedef Box2P2CJacobian<Grid, Scalar> LocalJacobian;
    typedef LeafP1TwoPhaseModel<Grid, Scalar, ProblemType, LocalJacobian> ThisLeafP1TwoPhaseModel;
    typedef Box2P2C<Grid, Scalar, VtkMultiWriter> ThisType;

    typedef typename ThisLeafP1TwoPhaseModel::FunctionType FunctionType;
    typedef typename ThisLeafP1TwoPhaseModel::FunctionType::RepresentationType VectorType;
    typedef typename ThisLeafP1TwoPhaseModel::OperatorAssembler::RepresentationType MatrixType;

    typedef typename Grid::LeafGridView GV;
    typedef typename GV::IndexSet IS;

    enum{m = 2};

    //////////////////////
    // Stuff required for the new newton method

    //! The traits class for the new newton method.
    struct NewtonTraits {
        typedef ScalarT                                                 Scalar;
        typedef typename ThisLeafP1TwoPhaseModel::FunctionType      Function;
        typedef typename ThisType::LocalJacobian                    LocalJacobian;
        typedef typename ThisLeafP1TwoPhaseModel::OperatorAssembler JacobianAssembler;
    };

    // HACK: traits for the domain of the problem. this is incomplete...
    struct DomainTraits {
        typedef ScalarT   Scalar;
        typedef GridT     Grid;
    };

    typedef NewNewtonMethod<ThisType> NewtonMethod;
    typedef TwoPTwoCNewtonController<NewtonMethod> NewtonController;

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

    Box2P2C(const Grid& grid, ProblemType& prob)
        : ThisLeafP1TwoPhaseModel(grid, prob)// (this->size) vectors
    { }

    void initial()
    {
        typedef typename Grid::Traits::template Codim<0>::Entity Entity;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename Grid::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

        enum {dim = Grid::dimension};
        enum {dimworld = Grid::dimensionworld};

        const GV& gridview(this->grid().leafView());

        this->localJacobian().outPressureN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMassFracAir = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMassFracWater = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        //        this->localJacobian().outPermeability = vtkMultiWriter->template createField<Scalar, 1>(this->size);


        this->localJacobian().clearVisited();
        // initialize switch flags
        this->localJacobian().setSwitchedLocal();
        this->localJacobian().setSwitchedLocalToGlobal();
        this->localJacobian().resetSwitchedLocal();

        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it)
        {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Entity& entity = *it;

            this->localJacobian().fvGeom.update(entity);

            const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<Scalar, Scalar, dim>::general(gt,1);

            int size = sfs.size();

            for (int i = 0; i < size; i++)
            {
                // get cell center in reference element
                const Dune::FieldVector<Scalar,dim>&local = sfs[i].position();

                // get global coordinate of cell center
                Dune::FieldVector<Scalar,dimworld> global = it->geometry().global(local);

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
            this->localJacobian().setLocalSolution(entity);
            this->localJacobian().updateStaticData(entity, this->localJacobian().u);
        }
        this->localJacobian().resetSwitchedLocal();
        this->localJacobian().resetSwitched();

        // set Dirichlet boundary conditions
        for (Iterator it = gridview.template begin<0>();
             it != eendit; ++it)
        {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Entity& entity = *it;

            const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<Scalar, Scalar, dim>::general(gt,1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().template assembleBC<LeafTag>(entity);

            IntersectionIterator
                endit = IntersectionIteratorGetter<Grid, LeafTag>::end(entity);

            for (IntersectionIterator is = IntersectionIteratorGetter<Grid,
                     LeafTag>::begin(entity); is!=endit; ++is)
                if (is->boundary())
                {
                    for (int i = 0; i < size; i++)
                        // handle subentities of this face
                        for (int j = 0; j < ReferenceElements<Scalar,dim>::general(gt).size(is->numberInSelf(), 1, sfs[i].codim()); j++)
                            if (sfs[i].entity()
                                == ReferenceElements<Scalar,dim>::general(gt).subEntity(is->numberInSelf(), 1,
                                                                                        j, sfs[i].codim()))
                            {
                                for (int equationNumber = 0; equationNumber<m; equationNumber++)
                                {
                                    if (this->localJacobian().bc(i)[equationNumber]
                                        == BoundaryConditions::dirichlet)
                                    {
                                        // get cell center in reference element
                                        Dune::FieldVector<Scalar,dim>
                                            local = sfs[i].position();

                                        // get global coordinate of cell center
                                        Dune::FieldVector<Scalar,dimworld>
                                            global = it->geometry().global(local);

                                        int globalId = this->vertexmapper.template map<dim>(entity, sfs[i].entity());
                                        FieldVector<int,m> dirichletIndex;
                                        FieldVector<BoundaryConditions::Flags, m>
                                            bctype = this->problem.bctype(global, entity, is, local);
                                        this->problem.dirichletIndex(global, entity, is,
                                                                     local, dirichletIndex);

                                        if (bctype[equationNumber] == BoundaryConditions::dirichlet)
                                        {
                                            FieldVector<Scalar,m>
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
        this->localJacobian().outPressureN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMassFracAir = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMassFracWater = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        //        this->localJacobian().outPermeability = vtkMultiWriter->template createField<Scalar, 1>(this->size);

        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);

        // execute newton method
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename Grid::Traits::template Codim<0>::Entity Entity;
        enum {dim = Grid::dimension};

        bool newtonLoop = false;
        while(!newtonLoop)
        {
            nextDt = this->localJacobian().getDt();
            NewtonMethod newton(*this); // *this means object itself (box2p2c)
            NewtonController newtonCtl(1e-7, 6);
            newtonLoop = newton.execute(*this, newtonCtl);
            nextDt = newtonCtl.suggestTimeStepSize(nextDt);
            this->localJacobian().setDt(nextDt);
            if(!newtonLoop){
                *this->u = *this->uOldTimeStep;
                this->localJacobian().resetPhaseState();
            }
            std::cout<<"timeStep resized to: "<<nextDt<<std::endl;
        }

        //        double Flux(0), Mass(0);
        //        Flux = this->computeFlux();
        //        Mass = totalCO2Mass();

        this->localJacobian().updatePhaseState(); // update variable oldPhaseState
        this->localJacobian().clearVisited();
        clearSwitched();
        //        updateState();                            // phase switch after each timestep
        //        std::cout << Flux << ", "<< Mass;

        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    virtual double totalCO2Mass() {
        typedef typename Grid::Traits::template Codim<0>::Entity Entity;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        enum {dim = Grid::dimension};
        enum {dimworld = Grid::dimensionworld};
        enum {gasPhase = 0, waterPhase = 1, bothPhases = 2};    // Phase state
        const GV& gridview(this->grid_.leafView());
        double totalMass = 0;
        double minSat = 1e100;
        double maxSat = -1e100;
        double minP  = 1e100;
        double maxP = -1e100;
        double minTe = 1e100;
        double maxTe = -1e100;
        double minX = 1e100;
        double maxX = -1e100;

        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it) {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Entity& entity = *it;

            FVElementGeometry<Grid> fvGeom;
            fvGeom.update(entity);

            const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<Scalar, Scalar, dim>::general(gt,
                                                                                1);
            int size = sfs.size();

            for (int i = 0; i < size; i++) {
                // get cell center in reference element
                const Dune::FieldVector<Scalar,dim>&local = sfs[i].position();

                // get global coordinate of cell center
                Dune::FieldVector<Scalar,dimworld> global = it->geometry().global(local);

                int globalId = this->vertexmapper.template map<dim>(entity,
                                                                    sfs[i].entity());

                int state;
                state = this->localJacobian().sNDat[globalId].phaseState;
                Scalar vol = fvGeom.subContVol[i].volume;
                Scalar poro = this->problem.soil().porosity(global, entity, local);

                Scalar rhoN = (*(this->localJacobian().outDensityN))[globalId];
                Scalar rhoW = (*(this->localJacobian().outDensityW))[globalId];
                Scalar satN = (*(this->localJacobian().outSaturationN))[globalId];
                Scalar satW = (*(this->localJacobian().outSaturationW))[globalId];
                Scalar xAW = (*(this->localJacobian().outMassFracAir))[globalId];
                Scalar xWN = (*(this->localJacobian().outMassFracWater))[globalId];
                Scalar xAN = 1 - xWN;
                //                Scalar pW = (*(this->u))[globalId][0];
                Scalar Te = (*(this->u))[globalId][2];
                Scalar mass = vol * poro * (satN * rhoN * xAN + satW * rhoW * xAW);



                minSat = std::min(minSat, satN);
                maxSat = std::max(maxSat, satN);
                //                minP = std::min(minP, pW);
                //                maxP = std::max(maxP, pW);
                minX = std::min(minX, xAW);
                maxX = std::max(maxX, xAW);
                minTe = std::min(minTe, Te);
                maxTe = std::max(maxTe, Te);

                totalMass += mass;
            }

        }

        // print minimum and maximum values
        std::cout << "nonwetting phase saturation: min = "<< minSat
                  << ", max = "<< maxSat << std::endl;
        std::cout << "wetting phase pressure: min = "<< minP
                  << ", max = "<< maxP << std::endl;
        std::cout << "mass fraction CO2: min = "<< minX
                  << ", max = "<< maxX << std::endl;
        std::cout << "mass: "<< totalMass << std::endl;

        return totalMass;
    }


    void restart(int restartNum)
    {
        typedef typename Grid::Traits::template Codim<0>::Entity Entity;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

        enum {dim = Grid::dimension};
        enum {dimworld = Grid::dimensionworld};

        const GV& gridview(this->grid().leafView());

        this->localJacobian().outPressureN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMassFracAir = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMassFracWater = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<Scalar, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<Scalar, 1>(this->size);

        int size = this->vertexmapper.size();
        Dune::BlockVector<FieldVector<double, m+1> > data(size);
        data=0;

        //        importFromDGF<GV>(data, "data", false);

        for (int i=0;i<size;i++)
        {
            for (int j=0;j<m;j++)
            {
                (*(this->u))[i][j]=data[i][j];
            }
        }

        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>();
             it != eendit; ++it)
        {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Entity& entity = *it;

            this->localJacobian().fvGeom.update(entity);

            const typename Dune::LagrangeShapeFunctionSetContainer<Scalar,Scalar,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<Scalar, Scalar, dim>::general(gt,1);
            int size = sfs.size();

            for (int i = 0; i < size; i++)
            {
                // get cell center in reference element
                const Dune::FieldVector<Scalar,dim>&local = sfs[i].position();

                // get global coordinate of cell center
                Dune::FieldVector<Scalar,dimworld> global = it->geometry().global(local);

                int globalId = this->vertexmapper.template map<dim>(entity,
                                                                    sfs[i].entity());

                // initialize cell concentration

                // initialize variable phaseState
                this->localJacobian().sNDat[globalId].phaseState = data[globalId][m];
                // initialize variable oldPhaseState
                this->localJacobian().sNDat[globalId].oldPhaseState = data[globalId][m];

            }
            this->localJacobian().clearVisited();
            this->localJacobian().setLocalSolution(entity);
            this->localJacobian().updateStaticData(entity, this->localJacobian().u);
        }
    }

    virtual void globalDefect(FunctionType& defectGlobal)
    {
        ThisLeafP1TwoPhaseModel::globalDefect(defectGlobal);
    }

    void solve()
    {
        return;
    }

    //    {
    //        Operator op(*(this->A));  // make operator out of matrix
    //        double red=1E-11;
    //
    //#ifdef HAVE_PARDISO
    ////    SeqPardiso<MatrixType,VectorType,VectorType> ilu0(*(this->A));
    //        pardiso.factorize(*(this->A));
    //        BiCGSTABSolver<VectorType> solver(op,pardiso,red,100,2);         // an inverse operator
    //    //    SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
    //        //LoopSolver<VectorType> solver(op, ilu0, red, 10, 2);
    //#else
    //        SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
    //
    //        //SeqIdentity<MatrixType,VectorType,VectorType> ilu0(*(this->A));// a precondtioner
    //        BiCGSTABSolver<VectorType> solver(op,ilu0,red,1000,1);         // an inverse operator
    //        //iteration numbers were 1e4 before
    //#endif
    //        InverseOperatorResult r;
    //        solver.apply(*(this->u), *(this->f), r);
    //
    //        return;
    //    }

    void updateState()
    {
        typedef typename Grid::Traits::template Codim<0>::Entity Entity;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename IntersectionIteratorGetter<Grid,LeafTag>::IntersectionIterator IntersectionIterator;

        enum {dim = Grid::dimension};
        enum {dimworld = Grid::dimensionworld};

        const GV& gridview(this->grid_.leafView());
        // iterate through leaf grid and evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>();
             it != eendit; ++it) {

            const Entity& entity = *it;
            this->localJacobian().fvGeom.update(entity);
            this->localJacobian().setLocalSolution(entity);
            this->localJacobian().computeElementData(entity);
            this->localJacobian().updateStaticData(entity, this->localJacobian().u);
        }

        return;
    }

    void clearSwitched()
    {
        this->localJacobian().resetSwitchedLocal();
        this->localJacobian().resetSwitched();
    }

    void setSwitchedLocalToGlobal()
    {
        this->localJacobian().setSwitchedLocalToGlobal();
    }

    bool checkSwitched()
    {
        return this->localJacobian().checkSwitched();
    }

    template<class MultiWriter>
    void addvtkfields (MultiWriter& writer)
    {
        //        BlockVector<FieldVector<Scalar, 1> > &xWN = *writer.template createField<Scalar, 1>(this->size);
        //        BlockVector<FieldVector<Scalar, 1> > &xAW = *writer.template createField<Scalar, 1>(this->size);
        //        BlockVector<FieldVector<Scalar, 1> > &satW = *writer.template createField<Scalar, 1>(this->size);

        //        writer.addScalarVertexFunction("nonwetting phase saturation", this->u, 1);
        writer.addScalarVertexFunction("pressure wetting phase", this->u, 0);

        //        writer.addVertexData(&satW,"wetting phase saturation");
        writer.addVertexData(this->localJacobian().outPressureN,"pressure non-wetting phase");
        writer.addVertexData(this->localJacobian().outCapillaryP,"capillary pressure");
        writer.addVertexData(this->localJacobian().outSaturationW,"saturation wetting phase");
        writer.addVertexData(this->localJacobian().outSaturationN,"saturation non-wetting phase");
        writer.addVertexData(this->localJacobian().outMassFracAir,"massfraction air in wetting phase");
        writer.addVertexData(this->localJacobian().outMassFracWater,"massfraction water in non-wetting phase");
        writer.addVertexData(this->localJacobian().outDensityW,"density wetting phase");
        writer.addVertexData(this->localJacobian().outDensityN,"density non-wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityW,"mobility wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityN,"mobility non-wetting phase");
        writer.addVertexData(this->localJacobian().outPhaseState,"phase state");
        //        writer.addVertexData(this->localJacobian().outPermeability,"permeability");
        //        writer.addVertexData(&xWN, "water in air");
        //        writer.addVertexData(&xAW, "dissolved air");
    }


    void vtkout (const char* name, int k)
    {

    }

    //    void vtkout (const char* name, int k)
    //    {
    //        VTKWriter<typename Grid::LeafGridView> vtkwriter(this->leafView());
    //        char fname[128];
    //        sprintf(fname,"%s-%05d",name,k);
    //      BlockVector<FieldVector<Scalar, 1> > xWN(this->size);
    //      BlockVector<FieldVector<Scalar, 1> > xAW(this->size);
    //        for (int i = 0; i < this->size; i++) {
    //            this->pW[i] = (*(this->u))[i][0];
    //            this->satN[i] = (*(this->u))[i][1];
    //            this->satW[i] = 1 - this->satN[i];
    //            //const FieldVector<Scalar, 4> parameters(this->problem.materialLawParameters
    //            //              (this->fvGeom.elementGlobal, e, this->fvGeom.elementLocal));
    //            //            this->pC[i] = this->problem.materialLaw().pC(this->satW[i], parameters);
    //            xWN[i] = this->problem.multicomp().xWN(this->pW[i], 283.15); //Achtung!! pW instead of pN!!!
    //            xAW[i] = this->problem.multicomp().xAW(this->pW[i], 283.15); //Achtung!! pW instead of pN!!!
    //        }
    //        vtkwriter.addVertexData(this->pW,"wetting phase pressure");
    //        vtkwriter.addVertexData(this->satW,"wetting phase saturation");
    //        vtkwriter.addVertexData(this->satN,"nonwetting phase saturation");
    //        vtkwriter.addVertexData(xWN, "water in air");
    //        vtkwriter.addVertexData(xAW, "dissolved air");
    //        vtkwriter.write(fname, VTKOptions::ascii);
    //    }
    //

    void setVtkMultiWriter(VtkMultiWriter *writer)
    { vtkMultiWriter = writer; }

protected:
    VtkMultiWriter *vtkMultiWriter;
    bool switchFlag;
};
}
#endif
