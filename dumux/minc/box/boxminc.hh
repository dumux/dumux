#ifndef DUNE_BOXMINC_HH
#define DUNE_BOXMINC_HH

#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>

//#include <dune/common/array.hh>      // defines simple array class
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
#include <dune/disc/functions/p0function.hh>
//#include <dune/disc/functions/p1function.hh>
#include "dumux/operators/p1operatorextended.hh"
#include <dune/disc/operators/boundaryconditions.hh>
/* #include <dune/disc/groundwater/groundwater.hh>
   #include <dune/disc/groundwater/p1groundwater.hh>
   #include <dune/disc/groundwater/p1groundwaterestimator.hh>
*/
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/istl/paamg/amg.hh>
#include "dumux/pardiso/pardiso.hh"
#include "dumux/pardiso/identity.hh"
#include "dumux/nonlinear/newtonmethod.hh"
#include "dumux/minc/mincmodel.hh"
#include "dumux/minc/mincproblem.hh"
#include "dumux/minc/box/boxmincjacobian.hh"

namespace Dune {

/**
   \brief Two phase model with Pw and Sn as primary unknowns

   This implements a two phase model with Pw and Sn as primary unknowns.
*/
template<class GridT, class Scalar, int numEq>
class BoxMinc :
        public LeafP1MincModel<GridT, Scalar, MincProblem<GridT, Scalar, numEq>,
                               BoxMincJacobian<GridT, Scalar, numEq>, numEq > {

public:

    //    typedef ScalarT Scalar;
    typedef GridT Grid;

    // define the problem type (also change the template argument above)
    typedef MincProblem<Grid, Scalar, numEq> ProblemType;

    // define the local Jacobian (also change the template argument above)
    typedef BoxMincJacobian<Grid, Scalar, numEq> LocalJacobian;

    typedef Dune::LeafP1MincModel<Grid, Scalar, ProblemType, LocalJacobian, numEq>
    LeafP1MincModel;

    typedef typename LeafP1MincModel::FunctionType FunctionType;

    typedef typename Grid::LeafGridView GV;
    typedef typename GV::IndexSet IS;

    //    enum {m = 4};

    typedef BoxMinc<GridT, Scalar, numEq> ThisType;
    typedef typename LeafP1MincModel::FunctionType::RepresentationType
    VectorType;
    typedef typename LeafP1MincModel::OperatorAssembler::RepresentationType
    MatrixType;
    typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
#ifdef HAVE_PARDISO
    SeqPardiso<MatrixType,VectorType,VectorType> pardiso;
#endif

    BoxMinc(const GridT& grid, ProblemType& prob) :
        LeafP1MincModel(grid, prob) {
    }

    void solve() {

        Operator op(*(this->A)); // make operator out of matrix
        double red=1E-8;

#ifdef HAVE_PARDISO
        //    SeqPardiso<MatrixType,VectorType,VectorType> ilu0(*(this->A));
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

        return;
    }

    void update(double& dt) {
        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);
        NewtonMethod<GridT, ThisType> newtonMethod(this->grid(), *this);
        newtonMethod.execute();
        dt = this->localJacobian().getDt();
        double upperMass, oldUpperMass;
        double totalMass = this->injected(upperMass, oldUpperMass);
        std::cout << totalMass << "\t"<< upperMass<< "\t"<< oldUpperMass
                  << "\t# totalMass, upperMass, oldUpperMass"<< std::endl;

        *(this->uOldTimeStep) = *(this->u);

        if (this->problem.exsolution)
            this->problem.updateExSol(dt, *(this->u));

        return;
    }

    virtual void globalDefect(FunctionType& defectGlobal) {
        typedef typename Grid::Traits::template Codim<0>::Entity Entity;
        typedef typename Grid::ctype DT;
        typedef typename GV::template Codim<0>::Iterator Iterator;
        enum {dim = Grid::dimension};
        typedef array<BoundaryConditions::Flags, numEq> BCBlockType;

        const GV& gridview(this->grid().leafView());
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
            bool old = true;
            this->localJacobian().updateVariableData(entity, this->localJacobian().uold, old);
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
                        essential[globalId].assign(BoundaryConditions::dirichlet);
                }
            }
        }

        for (typename std::vector<BCBlockType>::size_type i=0; i
                 <essential.size(); i++)
            for (int equationnumber = 0; equationnumber < numEq; equationnumber++) {
                if (essential[i][equationnumber] == BoundaryConditions::dirichlet)
                    (*defectGlobal)[i][equationnumber] = 0;
            }
    }
};

}
#endif
