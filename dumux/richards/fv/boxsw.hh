// $Id: boxpwsn.hh 726 2008-10-22 16:42:22Z bernd $

#ifndef DUNE_BOXPWSN_HH
#define DUNE_BOXPWSN_HH

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
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/paamg/amg.hh>
#include "dumux/pardiso/pardiso.hh"
#include "dumux/pardiso/identity.hh"
#include "dumux/nonlinear/newtonmethod.hh"
#include "dumux/richards/richardsmodel.hh"
#include "dumux/richards/richardsproblem.hh"
#include "dumux/richards/fv/boxswjacobian.hh"

namespace Dune {

/**
   \brief Two phase model with Pw and Sn as primary unknowns

   This implements a two phase model with Pw and Sn as primary unknowns.
*/
template<class G, class RT>
class BoxSw
    : public LeafP1TwoPhaseModel<G, RT, RichardsProblem<G, RT>, BoxSwJacobian<G, RT> >
{

public:
    // define the problem type (also change the template argument above)
    typedef RichardsProblem<G, RT> ProblemType;

    // define the local Jacobian (also change the template argument above)
    typedef BoxSwJacobian<G, RT> LocalJacobian;

    typedef LeafP1TwoPhaseModel<G, RT, ProblemType, LocalJacobian>
    ThisLeafP1TwoPhaseModel;

    typedef typename ThisLeafP1TwoPhaseModel::FunctionType FunctionType;

    typedef typename G::Traits::LeafIndexSet IS;

    enum {m = 1};

    typedef BoxSw<G, RT> ThisType;
    typedef typename ThisLeafP1TwoPhaseModel::FunctionType::RepresentationType VectorType;
    typedef typename ThisLeafP1TwoPhaseModel::OperatorAssembler::RepresentationType MatrixType;
    typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#if HAVE_MPI
#else
#ifdef HAVE_PARDISO
    SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
#endif
#endif

    BoxSw(const G& g, ProblemType& prob)
        : ThisLeafP1TwoPhaseModel(g, prob)
    {}

    virtual void solve()
    {
        Operator op(*(this->A));  // make operator out of matrix
        double red=1E-14;

#if HAVE_MPI
        // set up parallel solvers
        typedef typename G::Traits::GlobalIdSet::IdType GlobalIdType;
        typedef OwnerOverlapCopyExtendedCommunication<GlobalIdType,int> CommunicationType;
        Dune::IndexInfoFromGrid<GlobalIdType,int> indexinfo;
        (this->u).fillIndexInfoFromGrid(indexinfo);
        typedef Dune::OwnerOverlapCopyExtendedCommunication<GlobalIdType,int> CommunicationType;
        CommunicationType oocc(indexinfo,(this->grid).comm());
        int verbose=0;
        if ((this->grid).comm().rank() == 0)
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
        BiCGSTABSolver<VectorType> solver(op,pardiso,red,100,2);         // an inverse operator
#else
        SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
        BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,1);         // an inverse operator
#endif
        InverseOperatorResult r;
        solver.apply(*(this->u), *(this->f), r);
#endif

        return;
    }

    void update(double& dt) {
        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);
        NewtonMethod<G, ThisType> newtonMethod(this->grid(), *this);
        newtonMethod.execute();
        dt = this->localJacobian().getDt();
        double upperMass, oldUpperMass;
        double totalMass = this->injected(upperMass, oldUpperMass);
        std::cout << totalMass << "\t"<< upperMass<< "\t"<< oldUpperMass
                  << "\t"; //# totalMass, upperMass, oldUpperMass"<< std::endl;
        *(this->uOldTimeStep) = *(this->u);

        if (this->problem.exsolution)
            this->problem.updateExSol(dt, *(this->u));

    }

    virtual void vtkout(const char* name, int k) {
        VTKWriter<typename G::LeafGridView> vtkwriter(this->grid_.leafView());
        int size=this->vertexmapper.size();
        enum {dim = G::dimension};
        char fname[128];
        sprintf(fname, "%s-%05d", name, k);
        double minSat = 1e100;
        double maxSat = -1e100;
        if (this->problem.exsolution) {
            this->satEx.resize(size);
            this->satError.resize(size);
        }
        for (int i = 0; i < size; i++) {
            this->satW[i] = (*(this->u))[i][0];
            double satWI = this->satW[i];
            minSat = std::min(minSat, satWI);
            maxSat = std::max(maxSat, satWI);
            FieldVector<RT, dim> dummy (0);
            this->pC[i]  = this->problem.materialLaw().pC(this->satW[i], dummy, *(this->grid_.template leafbegin<0>()), dummy); // functions as long as we have constant properties for pc
            this->pW[i] = -this->pC[i]; // pw formulation
            if (this->problem.exsolution) {
                this->satEx[i]=this->problem.uExOutVertex(i, 1);
                this->satError[i]=this->problem.uExOutVertex(i, 2);
            }
        }
        vtkwriter.addVertexData(this->pW, "wetting phase pressure");
        vtkwriter.addVertexData(this->pC, "capillary pressure");
        vtkwriter.addVertexData(this->satW, "wetting phase saturation");
        if (this->problem.exsolution) {
            vtkwriter.addVertexData(this->satEx, "saturation, exact solution");
            vtkwriter.addVertexData(this->satError, "saturation error");
        }
        vtkwriter.write(fname, VTKOptions::ascii);
        std::cout << "wetting phase saturation: min = "<< minSat
                  << ", max = "<< maxSat << std::endl;
        if (minSat< -0.5 || maxSat > 1.5)DUNE_THROW(MathError, "Saturation exceeds range.");
    }
};

}
#endif
