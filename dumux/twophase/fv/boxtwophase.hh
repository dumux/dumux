// $Id$

#ifndef DUNE_BOXTWOPHASE_HH
#define DUNE_BOXTWOPHASE_HH

#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>

#include "dumux/nonlinear/newtonmethod.hh"
#include "dumux/twophase/twophasemodel_deprecated.hh"
#include "dumux/twophase/twophaseproblem_deprecated.hh"
#include "dumux/twophase/fv/boxtwophasejacobian.hh"

namespace Dune
{
/** \todo Please doc me! */

template<class Grid, class Scalar>
class BoxTwoPhase
    : public LeafP1TwoPhaseModel<Grid, Scalar, TwoPhaseProblem<Grid, Scalar>,
                                 BoxTwoPhaseLocalJacobian<Grid, Scalar> >
{
public:
    // define the problem type (also change the template argument above)
    typedef TwoPhaseProblem<Grid, Scalar> ProblemType;

    // define the local Jacobian (also change the template argument above)
    typedef BoxTwoPhaseLocalJacobian<Grid, Scalar> LocalJacobian;

    typedef LeafP1TwoPhaseModel<Grid, Scalar, ProblemType, LocalJacobian> LeafP1TwoPhaseModel;

    typedef BoxTwoPhase<Grid, Scalar> ThisType;

    BoxTwoPhase(const Grid& grid, ProblemType& prob)
        : LeafP1TwoPhaseModel(grid, prob)
    {     }

    void solve()
    {
        typedef typename LeafP1TwoPhaseModel::FunctionType::RepresentationType VectorType;
        typedef typename LeafP1TwoPhaseModel::OperatorAssembler::RepresentationType MatrixType;
        typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator;

        Operator op(*(this->A));  // make operator out of matrix
        double red=1E-14;
        SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
        BiCGSTABSolver<VectorType> solver(op,ilu0,red,10000,0);         // an inverse operator
        InverseOperatorResult r;
        solver.apply(*(this->u), *(this->f), r);

        return;
    }

    void update (double& dt)
    {
        this->localJacobian().setDt(dt);
        this->localJacobian().setOldSolution(this->uOldTimeStep);
        NewtonMethod<Grid, ThisType> newtonMethod(this->grid, *this);
        newtonMethod.execute();
        *(this->uOldTimeStep) = *(this->u);

        return;
    }

};

}
#endif
