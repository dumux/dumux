// $Id$

#ifndef DUNE_BOX1P2C_HH
#define DUNE_BOX1P2C_HH

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
#include "dumux/pardiso/identity.hh"
#include "dumux/nonlinear/newtonmethod.hh"
#include "dumux/1p2c/onephasemodel.hh"
#include "dumux/1p2c/1p2cproblem.hh"
#include "dumux/1p2c/fv/box1p2cjacobian.hh"

#include "dumux/nonlinear/new_newtonmethod.hh"
#include "dumux/nonlinear/new_newtoncontroller.hh"

namespace Dune
{
/**
   \brief Isothermal one phase two component model with P and X as primary unknowns

   This implements an isothermal one phase two component model with P and X as primary unknowns
*/
template<class GridT, class ScalarRT, class VtkMultiWriter>
class Box1P2C
    : public LeafP1OnePhaseModel<GridT, ScalarRT, OnePTwoCProblem<GridT, ScalarRT>, Box1P2CJacobian<GridT, ScalarRT> >
{
public:
    // define the problem type (also change the template argument above)
    typedef OnePTwoCProblem<GridT, ScalarRT> ProblemType;

    // define the local Jacobian (also change the template argument above)
    typedef Box1P2CJacobian<GridT, ScalarRT> LocalJacobian;
    typedef LeafP1OnePhaseModel<GridT, ScalarRT, ProblemType, LocalJacobian> ThisLeafP1OnePhaseModel;
    typedef Box1P2C<GridT, ScalarRT, VtkMultiWriter> ThisType;

    typedef typename ThisLeafP1OnePhaseModel::FunctionType FunctionType;
    typedef typename ThisLeafP1OnePhaseModel::FunctionType::RepresentationType VectorType;
    typedef typename ThisLeafP1OnePhaseModel::OperatorAssembler::RepresentationType MatrixType;
    typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#ifdef HAVE_PARDISO
    SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
#endif
    typedef typename GridT::LeafGridView GV;


    enum{numEq = 2};

    //////////////////////
    // Stuff required for the new newton method

    //! The traits class for the new newton method.
    struct NewtonTraits {
        typedef ScalarRT                                            Scalar;
        typedef typename ThisLeafP1OnePhaseModel::FunctionType      Function;
        typedef typename ThisType::LocalJacobian                    LocalJacobian;
        typedef typename ThisLeafP1OnePhaseModel::OperatorAssembler JacobianAssembler;
    };

    // HACK: traits for the domain of the problem. this is incomplete...
    struct DomainTraits {
        typedef ScalarRT   Scalar;
        typedef GridT    Grid;
    };

    typedef NewNewtonMethod<ThisType> NewtonMethod;
    typedef Dune::NewtonController<NewtonMethod> NewtonController;

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

    Box1P2C(const GridT& g, ProblemType& prob)
        : ThisLeafP1OnePhaseModel(g, prob)// (this->size) vectors
    { }

    void initial()
    {
        typedef typename GridT::Traits::template Codim<0>::Entity Element;
        typedef typename GridT::ctype CoordScalar;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename GridT::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

        enum {dim = GridT::dimension};
        enum {dimworld = GridT::dimensionworld};

        const GV& gridview(this->grid().leafView());


        this->localJacobian().clearVisited();


        // iterate through leaf grid an evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>(); it
                 != eendit; ++it)
        {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get entity
            const Element& element = *it;

            this->localJacobian().fvGeom.update(element);

            const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,ScalarRT,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<CoordScalar, ScalarRT, dim>::general(gt,1);

            int size = sfs.size();

            for (int i = 0; i < size; i++)
            {
                // get cell center in reference element
                const Dune::FieldVector<CoordScalar,dim>&local = sfs[i].position();

                // get global coordinate of cell center
                Dune::FieldVector<CoordScalar,dimworld> global = it->geometry().global(local);

                int globalId = this->vertexmapper.template map<dim>(element,
                                                                    sfs[i].entity());

                // initialize cell concentration
                (*(this->u))[globalId] = this->problem.initial(
                                                               global, element, local);
            }
            this->localJacobian().clearVisited();
            this->localJacobian().setLocalSolution(element);
            this->localJacobian().updateStaticData(element, this->localJacobian().u);
        }

        // set Dirichlet boundary conditions
        for (Iterator it = gridview.template begin<0>();
             it != eendit; ++it)
        {
            // get geometry type
            Dune::GeometryType gt = it->geometry().type();

            // get element
            const Element& element = *it;

            const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,ScalarRT,dim>::value_type
                &sfs=Dune::LagrangeShapeFunctions<CoordScalar, ScalarRT, dim>::general(gt,1);
            int size = sfs.size();

            // set type of boundary conditions
            this->localJacobian().assembleBoundaryCondition(element);

            typedef typename GridT::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

            IntersectionIterator endit = element.ileafend();
            for (IntersectionIterator is = element.ileafbegin(); is!=endit; ++is)
                if (is->boundary())
                {
                    for (int i = 0; i < size; i++)
                        // handle subentities of this face
                        for (int j = 0; j < ReferenceElements<CoordScalar,dim>::general(gt).size(is->indexInInside(), 1, sfs[i].codim()); j++)
                            if (sfs[i].entity()
                                == ReferenceElements<CoordScalar,dim>::general(gt).subEntity(is->indexInInside(), 1,
                                                                                             j, sfs[i].codim()))
                            {
                                for (int equationNumber = 0; equationNumber<numEq; equationNumber++)
                                {
                                    if (this->localJacobian().bc(i)[equationNumber]
                                        == BoundaryConditions::dirichlet)
                                    {
                                        // get cell center in reference element
                                        Dune::FieldVector<CoordScalar,dim>
                                            local = sfs[i].position();

                                        // get global coordinate of cell center
                                        Dune::FieldVector<CoordScalar,dimworld>
                                            global = it->geometry().global(local);

                                        int globalId = this->vertexmapper.template map<dim>(element, sfs[i].entity());
                                        FieldVector<int,numEq> dirichletIndex;
                                        FieldVector<BoundaryConditions::Flags, numEq>
                                            bctype = this->problem.bctype(global, element, is, local);
                                        this->problem.dirichletIndex(global, element, is,
                                                                     local, dirichletIndex);

                                        if (bctype[equationNumber] == BoundaryConditions::dirichlet)
                                        {
                                            FieldVector<ScalarRT,numEq>
                                                ghelp = this->problem.g(
                                                                        global, element, is,
                                                                        local);
                                            (*(this->u))[globalId][dirichletIndex[equationNumber]]
                                                = ghelp[dirichletIndex[equationNumber]];
                                        }
                                    }
                                }
                            }
                }
            this->localJacobian().setLocalSolution(element);
            for (int i = 0; i < size; i++)
                this->localJacobian().updateVariableData(element, this->localJacobian().u, i, false);

        }

        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    void updateModel (double& dt, double& nextDt)
    {

        this->localJacobian.setDt(dt);
        this->localJacobian.setOldSolution(this->uOldTimeStep);

        // execute newton method
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename GridT::Traits::template Codim<0>::Entity Element;
        typedef typename GridT::ctype CoordScalar;
        enum {dim = GridT::dimension};

        bool newtonLoop = false;
        while(!newtonLoop)
        {
            nextDt = this->localJacobian.getDt();
            NewtonMethod newton(*this); // *this means object itself (box1p2c)
            NewtonController newtonCtl(1e-7, 6);
            newtonLoop = newton.execute(*this, newtonCtl);
            nextDt = newtonCtl.suggestTimeStepSize(nextDt);
            this->localJacobian.setDt(nextDt);
            if(!newtonLoop){
                *this->u = *this->uOldTimeStep;
                this->localJacobian.resetPhaseState();
            }
            std::cout<<"timeStep resized to: "<<nextDt<<std::endl;
        }


        this->localJacobian.clearVisited();

        *(this->uOldTimeStep) = *(this->u);

        return;
    }

    //    void restart()
    //    {
    //        typedef typename GridT::Traits::template Codim<0>::Entity Element;
    //        typedef typename GridT::ctype CoordScalar;
    //        typedef typename GV::template Codim<0>::Iterator Iterator;
    //        typedef typename IntersectionIteratorGetter<GridT,LeafTag>::IntersectionIterator IntersectionIterator;
    //
    //        enum {dim = GridT::dimension};
    //        enum {dimworld = GridT::dimensionworld};
    //
    //        const GV& gridview(this->grid().leafView());
    //
    //        this->localJacobian.outPressureN = vtkMultiWriter->template createField<ScalarRT, 1>(this->size);
    //        this->localJacobian.outCapillaryP = vtkMultiWriter->template createField<ScalarRT, 1>(this->size);
    //        this->localJacobian.outSaturationW = vtkMultiWriter->template createField<ScalarRT, 1>(this->size);
    //        this->localJacobian.outSaturationN = vtkMultiWriter->template createField<ScalarRT, 1>(this->size);
    //        this->localJacobian.outMassFracAir = vtkMultiWriter->template createField<ScalarRT, 1>(this->size);
    //        this->localJacobian.outMassFracWater = vtkMultiWriter->template createField<ScalarRT, 1>(this->size);
    //        this->localJacobian.outDensityW = vtkMultiWriter->template createField<ScalarRT, 1>(this->size);
    //        this->localJacobian.outDensityN = vtkMultiWriter->template createField<ScalarRT, 1>(this->size);
    //        this->localJacobian.outMobilityW = vtkMultiWriter->template createField<ScalarRT, 1>(this->size);
    //        this->localJacobian.outMobilityN = vtkMultiWriter->template createField<ScalarRT, 1>(this->size);
    //        this->localJacobian.outPhaseState = vtkMultiWriter->template createField<ScalarRT, 1>(this->size);
    //
    //        int size = this->vertexmapper.size();
    //        Dune::BlockVector<FieldVector<double, numEq+1> > data(size);
    //        data=0;
    //
    ////        importFromDGF<GV>(data, "data", false);
    //
    //        for (int i=0;i<size;i++)
    //        {
    //            for (int j=0;j<numEq;j++)
    //            {
    //                (*(this->u))[i][j]=data[i][j];
    //            }
    //        }
    //
    //        // iterate through leaf grid an evaluate c0 at cell centere
    //        Iterator eendit = gridview.template end<0>();
    //        for (Iterator it = gridview.template begin<0>();
    //            it != eendit; ++it)
    //        {
    //            // get geometry type
    //            Dune::GeometryType gt = it->geometry().type();
    //
    //            // get element
    //            const Element& element = *it;
    //
    //            this->localJacobian.fvGeom.update(element);
    //
    //            const typename Dune::LagrangeShapeFunctionSetContainer<CoordScalar,ScalarRT,dim>::value_type
    //                    &sfs=Dune::LagrangeShapeFunctions<CoordScalar, ScalarRT, dim>::general(gt,1);
    //            int size = sfs.size();
    //
    //            for (int i = 0; i < size; i++)
    //            {
    //                // get cell center in reference element
    //                const Dune::FieldVector<CoordScalar,dim>&local = sfs[i].position();
    //
    //                // get global coordinate of cell center
    //                Dune::FieldVector<CoordScalar,dimworld> global = it->geometry().global(local);
    //
    //                int globalId = this->vertexmapper.template map<dim>(element,
    //                        sfs[i].entity());
    //
    //                // initialize cell concentration
    //
    //                // initialize variable phaseState
    //                this->localJacobian.sNDat[globalId].phaseState = data[globalId][numEq];
    //                // initialize variable oldPhaseState
    //                this->localJacobian.sNDat[globalId].oldPhaseState = data[globalId][numEq];
    //
    //            }
    //            this->localJacobian.clearVisited();
    //            this->localJacobian.setLocalSolution(element);
    //            this->localJacobian.updateStaticData(element, this->localJacobian.u);
    //        }
    //    }

    virtual void globalDefect(FunctionType& defectGlobal)
    {
        ThisLeafP1OnePhaseModel::globalDefect(defectGlobal);
    }

    void solve()
    {
        Operator op(*(this->A));  // make operator out of matrix
        double red=1E-14;

#ifdef HAVE_PARDISO
        pardiso.factorize(*(this->A));
        BiCGSTABSolver<VectorType> solver(op,pardiso,red,100,2);         // an inverse operator
#else
        SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
        BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,1);         // an inverse operator
#endif
        InverseOperatorResult r;
        solver.apply(*(this->u), *(this->f), r);

        return;
    }

    void updateState()
    {
        typedef typename GridT::Traits::template Codim<0>::Entity Element;
        typedef typename GridT::ctype CoordScalar;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        typedef typename GridT::template Codim<0>::LeafIntersectionIterator IntersectionIterator;

        enum {dim = GridT::dimension};
        enum {dimworld = GridT::dimensionworld};

        const GV& gridview(this->grid_.leafView());
        // iterate through leaf grid and evaluate c0 at cell center
        Iterator eendit = gridview.template end<0>();
        for (Iterator it = gridview.template begin<0>();
             it != eendit; ++it) {

            const Element& element = *it;
            this->localJacobian.fvGeom.update(element);
            this->localJacobian.setLocalSolution(element);
            this->localJacobian.computeElementData(element);
            this->localJacobian.updateStaticData(element, this->localJacobian.u);
        }

        return;
    }



    template<class MultiWriter>
    void addvtkfields (MultiWriter& writer)
    {

    }

    void vtkout (const char* name, int k)
    {
        VTKWriter<typename GridT::LeafGridView> vtkwriter(this->grid_.leafView());
        char fname[128];
        sprintf(fname,"%s-%05d",name,k);
        BlockVector<FieldVector<ScalarRT, 1> > p(this->size);
        BlockVector<FieldVector<ScalarRT, 1> > x(this->size);
        for (int i = 0; i < this->size; i++) {
            p[i] = (*(this->u))[i][0];
            x[i] = (*(this->u))[i][1];
            //const FieldVector<ScalarRT, 4> parameters(this->problem.materialLawParameters
            //              (this->fvGeom.cellGlobal, e, this->fvGeom.cellLocal));
            //            this->pC[i] = this->problem.materialLaw().pC(this->satW[i], parameters);
        }

        vtkwriter.addVertexData(p, "pressure");
        vtkwriter.addVertexData(x, "mole fraction");
        vtkwriter.write(fname, VTKOptions::ascii);
    }


    void setVtkMultiWriter(VtkMultiWriter *writer)
    { vtkMultiWriter = writer; }

protected:
    VtkMultiWriter *vtkMultiWriter;

};
}

#endif
