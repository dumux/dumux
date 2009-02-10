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
#include "dumux/operators/p1operatorextended.hh"
#include "dumux/operators/owneroverlapcopyextended.hh"
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

namespace Dune
{

/**
 \brief Two phase model with Pw and Sn as primary unknowns

 This implements a two phase model with Pw and Sn as primary unknowns.
 */
template<class ThisGrid, class ThisScalar, class VtkMultiWriter>
class BoxPwSnTe: public LeafP1TwoPhaseModel<ThisGrid, ThisScalar,
        TwoPhaseHeatProblem<ThisGrid, ThisScalar> , BoxPwSnTeJacobian<ThisGrid,
                ThisScalar> >
{

public:
    // define the problem type (also change the template argument above)
    typedef TwoPhaseHeatProblem<ThisGrid, ThisScalar> ProblemType;

    // define the local Jacobian (also change the template argument above)
    typedef BoxPwSnTeJacobian<ThisGrid, ThisScalar> LocalJacobian;
    typedef LeafP1TwoPhaseModel<ThisGrid, ThisScalar, ProblemType,
            LocalJacobian> ThisLeafP1TwoPhaseModel;
typedef    typename ThisLeafP1TwoPhaseModel::FunctionType FunctionType;

    typedef typename ThisGrid::LeafGridView GridView;
    typedef typename GridView::IndexSet IS;

    enum
    {   numEq = 3};

    typedef BoxPwSnTe<ThisGrid, ThisScalar, VtkMultiWriter> ThisType;
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
    struct NewtonTraits
    {
        typedef ThisScalar Scalar;
        typedef typename ThisLeafP1TwoPhaseModel::FunctionType Function;
        typedef typename ThisType::LocalJacobian LocalJacobian;
        typedef typename ThisLeafP1TwoPhaseModel::OperatorAssembler JacobianAssembler;
    };

    // HACK: traits for the domain of the problem. this is incomplete...
    struct DomainTraits
    {
        typedef ThisScalar Scalar;
        typedef ThisGrid Grid;
    };

    typedef NewNewtonMethod<ThisType> NewtonMethod;
    typedef TwoPTwoCNINewtonController<NewtonMethod> NewtonController;

    typedef typename NewtonTraits::Function Function;
    Function &currentSolution()
    {   return this->u;};

    LocalJacobian &getLocalJacobian()
    {   return this->localJacobian();}

    typedef typename NewtonTraits::JacobianAssembler JacobianAssembler;
    JacobianAssembler &jacobianAssembler()
    {   return this->A;}
    // End of stuff for new newton method
    //////////////////////

    BoxPwSnTe(const ThisGrid& grid, ProblemType& prob)
    : ThisLeafP1TwoPhaseModel(grid, prob)
    {}

    void initial()
    {
        typedef typename ThisGrid::Traits::template Codim<0>::Entity Element;
        typedef typename ThisGrid::ctype CoordScalar;
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        typedef typename IntersectionIteratorGetter<ThisGrid,LeafTag>::IntersectionIterator IntersectionIterator;

        enum
        {   dim = ThisGrid::dimension};
        enum
        {   dimworld = ThisGrid::dimensionworld};

        const GridView& gridview(this->_grid.leafView());

        this->localJacobian().outPressureN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outTemperature = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);

        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator eIt = gridview.template begin<0>(); eIt
                != eendit; ++eIt)
        {
            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get element
            const Element& element = *eIt;

            this->localJacobian().fvGeom.update(element);

            const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,ThisScalar,dim>::value_type
            &sfs=Dune::LagrangeShapeFunctions<CoordScalar, ThisScalar, dim>::general(gt,1);
            int size = sfs.size();

            for (int idx = 0; idx < size; idx++)
            {
                // get cell center in reference element
                const Dune::FieldVector<CoordScalar,dim>&localPos = sfs[idx].position();

                // get global coordinate of cell center
                Dune::FieldVector<CoordScalar,dimworld> globalPos = eIt->geometry().global(localPos);

                int globalId = this->vertexmapper.template map<dim>(element,
                        sfs[idx].entity());

                // initialize cell concentration
                (*(this->u))[globalId] = this->problem.initial(globalPos, element, localPos);
            }
            this->localJacobian().clearVisited();
            this->localJacobian().initiateStaticData(element);
        }

        // set Dirichlet boundary conditions
        for (Iterator eIt = gridview.template begin<0>(); eIt
                != eendit; ++eIt)
        {
            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get element
            const Element& element = *eIt;

            const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,ThisScalar,dim>::value_type
            &sfs=Dune::LagrangeShapeFunctions<CoordScalar, ThisScalar, dim>::general(gt,
                    1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().template assembleBC<LeafTag>(element);

            IntersectionIterator
            endit = IntersectionIteratorGetter<ThisGrid, LeafTag>::end(element);
            for (IntersectionIterator isIt = IntersectionIteratorGetter<ThisGrid,
                    LeafTag>::begin(element); isIt!=endit; ++isIt)
            if (isIt->boundary())
            {
                for (int idx = 0; idx < size; idx++)
                // handle subentities of this face
                for (int j = 0; j < ReferenceElements<CoordScalar,dim>::general(gt).size(isIt->numberInSelf(), 1, sfs[idx].codim()); j++)
                if (sfs[idx].entity()
                        == ReferenceElements<CoordScalar,dim>::general(gt).subEntity(isIt->numberInSelf(), 1,
                                j, sfs[idx].codim()))
                {
                    for (int equationNumber = 0; equationNumber<numEq; equationNumber++)
                    {
                        if (this->localJacobian().bc(idx)[equationNumber]
                                == BoundaryConditions::dirichlet)
                        {
                            // get cell center in reference element
                            Dune::FieldVector<CoordScalar,dim>
                            localPos = sfs[idx].position();

                            // get global coordinate of cell center
                            Dune::FieldVector<CoordScalar,dimworld>
                            globalPos = eIt->geometry().global(localPos);

                            int
                            globalId = this->vertexmapper.template map<dim>(
                                    element, sfs[idx].entity());
                            FieldVector<int,numEq> dirichletIndex;
                            FieldVector<BoundaryConditions::Flags, numEq>
                            bctype = this->problem.bctype(
                                    globalPos, element, isIt,
                                    localPos);
                            this->problem.dirichletIndex(globalPos, element, isIt,
                                    localPos, dirichletIndex);

                            if (bctype[equationNumber]
                                    == BoundaryConditions::dirichlet)
                            {
                                FieldVector<ThisScalar,numEq>
                                ghelp = this->problem.g(
                                        globalPos, element, isIt,
                                        localPos);
                                (*(this->u))[globalId][dirichletIndex[equationNumber]]
                                = ghelp[dirichletIndex[equationNumber]];
                            }
                        }
                    }
                }
            }
            this->localJacobian().setLocalSolution(element);
            for (int idx = 0; idx < size; idx++)
            this->localJacobian().updateVariableData(element, this->localJacobian().u, idx, false);

        }

        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    void restart(int restartNum=0)
    {
        typedef typename ThisGrid::Traits::template Codim<0>::Entity Element;
        typedef typename ThisGrid::ctype CoordScalar;
        typedef typename GridView::template Codim<0>::Iterator Iterator;
        typedef typename IntersectionIteratorGetter<ThisGrid,LeafTag>::IntersectionIterator IntersectionIterator;

        enum
        {   dim = ThisGrid::dimension};
        enum
        {   dimworld = ThisGrid::dimensionworld};

        const GridView& gridview(this->_grid.leafView());

        this->localJacobian().outPressureN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outCapillaryP = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outSaturationN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outTemperature = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outPhaseState = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);

        int size = this->vertexmapper.size();
        Dune::BlockVector<FieldVector<double, numEq> > data(size);
        data=0;

        std::string restartFileName;
        restartFileName = (boost::format("data-%05d")
                %restartNum).str();
        importFromDGF<GridView>(data, restartFileName, false);

        for (int globalIdx=0;globalIdx<size;globalIdx++)
        {
            for (int eq=0;eq<numEq;eq++)
            {
                (*(this->u))[globalIdx][eq]=data[globalIdx][eq];
            }
        }

        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator eIt = gridview.template begin<0>(); eIt
                != eendit; ++eIt)
        {
            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get element
            const Element& element = *eIt;

            this->localJacobian().fvGeom.update(element);

            const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,ThisScalar,dim>::value_type
            &sfs=Dune::LagrangeShapeFunctions<CoordScalar, ThisScalar, dim>::general(gt,1);
            int size = sfs.size();

            for (int idx = 0; idx < size; idx++)
            {
                // get cell center in reference element
                const Dune::FieldVector<CoordScalar,dim>&localPos = sfs[idx].position();

                // get global coordinate of cell center
                Dune::FieldVector<CoordScalar,dimworld> globalPos = eIt->geometry().global(localPos);

                int globalId = this->vertexmapper.template map<dim>(element,
                        sfs[idx].entity());
            }
            this->localJacobian().clearVisited();
            this->localJacobian().initiateStaticData(element);
        }

        // set Dirichlet boundary conditions
        for (Iterator eIt = gridview.template begin<0>(); eIt
                != eendit; ++eIt)
        {
            // get geometry type
            Dune::GeometryType gt = eIt->geometry().type();

            // get element
            const Element& element = *eIt;

            const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,ThisScalar,dim>::value_type
            &sfs=Dune::LagrangeShapeFunctions<CoordScalar, ThisScalar, dim>::general(gt,
                    1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().template assembleBC<LeafTag>(element);

            IntersectionIterator
            endit = IntersectionIteratorGetter<ThisGrid, LeafTag>::end(element);
            for (IntersectionIterator isIt = IntersectionIteratorGetter<ThisGrid,
                    LeafTag>::begin(element); isIt!=endit; ++isIt)
            if (isIt->boundary())
            {
                for (int idx = 0; idx < size; idx++)
                // handle subentities of this face
                for (int j = 0; j < ReferenceElements<CoordScalar,dim>::general(gt).size(isIt->numberInSelf(), 1, sfs[idx].codim()); j++)
                if (sfs[idx].entity()
                        == ReferenceElements<CoordScalar,dim>::general(gt).subEntity(isIt->numberInSelf(), 1,
                                j, sfs[idx].codim()))
                {
                    for (int equationNumber = 0; equationNumber<numEq; equationNumber++)
                    {
                        if (this->localJacobian().bc(idx)[equationNumber]
                                == BoundaryConditions::dirichlet)
                        {
                            // get cell center in reference element
                            Dune::FieldVector<CoordScalar,dim>
                            localPos = sfs[idx].position();

                            // get global coordinate of cell center
                            Dune::FieldVector<CoordScalar,dimworld>
                            globalPos = eIt->geometry().global(localPos);

                            int
                            globalId = this->vertexmapper.template map<dim>(
                                    element, sfs[idx].entity());
                            FieldVector<int,numEq> dirichletIndex;
                            FieldVector<BoundaryConditions::Flags, numEq>
                            bctype = this->problem.bctype(
                                    globalPos, element, isIt,
                                    localPos);
                            this->problem.dirichletIndex(globalPos, element, isIt,
                                    localPos, dirichletIndex);

                            if (bctype[equationNumber]
                                    == BoundaryConditions::dirichlet)
                            {
                                FieldVector<ThisScalar,numEq>
                                ghelp = this->problem.g(
                                        globalPos, element, isIt,
                                        localPos);
                                (*(this->u))[globalId][dirichletIndex[equationNumber]]
                                = ghelp[dirichletIndex[equationNumber]];
                            }
                        }
                    }
                }
            }
            this->localJacobian().setLocalSolution(element);
            for (int idx = 0; idx < size; idx++)
            this->localJacobian().updateVariableData(element, this->localJacobian().u, idx, false);

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
        Operator op(*(this->A)); // make operator out of matrix
        double red=1E-9;

#ifdef HAVE_PARDISO
        //    SeqPardiso<MatrixType,VectorType,VectorType> ilu0(*(this->A));
        pardiso.factorize(*(this->A));
        BiCGSTABSolver<VectorType> solver(op,pardiso,red,100,2); // an inverse operator
        //    SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
        //LoopSolver<VectorType> solver(op, ilu0, red, 10, 2);
#else
        SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner

        //SeqIdentity<MatrixType,VectorType,VectorType> ilu0(*(this->A));// a precondtioner
        BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,1); // an inverse operator
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
        this->localJacobian().outDensityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outDensityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityW = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        this->localJacobian().outMobilityN = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);
        //        this->localJacobian().outPermeability = vtkMultiWriter->template createField<ThisScalar, 1>(this->size);

        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);

        typedef typename GridView::template Codim<0>::Iterator Iterator;
        typedef typename ThisGrid::Traits::template Codim<0>::Entity Element;
        typedef typename ThisGrid::ctype CoordScalar;
        enum
        {   dim = ThisGrid::dimension};

        bool newtonLoop = false;
        while(!newtonLoop)
        {
            nextDt = this->localJacobian().getDt();
            NewtonMethod newton(*this);
            NewtonController newtonCtl;
            newtonLoop = newton.execute(*this, newtonCtl);
            nextDt = newtonCtl.suggestTimeStepSize(nextDt);
            this->localJacobian().setDt(nextDt);
            if(!newtonLoop)
            {
                *this->u = *this->uOldTimeStep;
            }
            std::cout<<"timeStep resized to: "<<nextDt<<std::endl;
        }

        this->localJacobian().clearVisited();
        //        std::cout << Flux << ", "<< Mass;

        *(this->uOldTimeStep) = *(this->u);

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
        writer.addVertexData(this->localJacobian().outPressureN,"pressure non-wetting phase");
        writer.addVertexData(this->localJacobian().outCapillaryP,"capillary pressure");
        writer.addVertexData(this->localJacobian().outTemperature,"temperature");
        writer.addVertexData(this->localJacobian().outSaturationW,"saturation wetting phase");
        writer.addVertexData(this->localJacobian().outSaturationN,"saturation non-wetting phase");
        writer.addVertexData(this->localJacobian().outDensityW,"density wetting phase");
        writer.addVertexData(this->localJacobian().outDensityN,"density non-wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityW,"mobility wetting phase");
        writer.addVertexData(this->localJacobian().outMobilityN,"mobility non-wetting phase");
    }

    void setVtkMultiWriter(VtkMultiWriter *writer)
    {   vtkMultiWriter = writer;}

protected:
    VtkMultiWriter *vtkMultiWriter;
};

}
#endif
