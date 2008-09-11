// $Id$ 

#ifndef DUNE_FENONLINEARPARABOLIC_HH
#define DUNE_FENONLINEARPARABOLIC_HH

#include<dune/istl/operators.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/preconditioners.hh>

#include "dumux/parabolic/nonlinearparabolic.hh"
#include "dumux/parabolic/parabolicproblem.hh"
#include "dumux/parabolic/fe/p1parabolicjacobian.hh"

namespace Dune
{
  template<class G, class RT>
  class FENonlinearParabolic 
  : public LeafP1NonlinearParabolic<G, RT, ParabolicProblem<G, RT>, 
                                 P1ParabolicLocalJacobian<G, RT>, 1>
  {
  public:
	// define the problem type:
	typedef ParabolicProblem<G, RT> ProblemType;
	
	// define the local Jacobian
	typedef P1ParabolicLocalJacobian<G, RT> LocalJacobian;
	
	typedef LeafP1NonlinearParabolic<G, RT, ProblemType, LocalJacobian, 1> LeafP1NonlinearParabolic;

	typedef FENonlinearParabolic<G, RT> ThisType;
	
	FENonlinearParabolic(const G& g, ProblemType& prob) 
	: LeafP1NonlinearParabolic(g, prob)
	{ 	}

	void solve() 
	{
		typedef typename LeafP1NonlinearParabolic::FunctionType::RepresentationType VectorType;
		typedef typename LeafP1NonlinearParabolic::OperatorAssembler::RepresentationType MatrixType;
		typedef MatrixAdapter<MatrixType,VectorType,VectorType> Operator; 
		
		Operator op(*(this->A));  // make operator out of matrix
		double red=1E-14;
		SeqILU0<MatrixType,VectorType,VectorType> ilu0(*(this->A),1.0);// a precondtioner
		CGSolver<VectorType> solver(op,ilu0,red,10000,0);         // an inverse operator 
		InverseOperatorResult r;
		solver.apply(*(this->u), *(this->f), r);
		
		return;		
	}

	void update (double& dt)
	{
		this->localJacobian.setDt(dt);
		this->localJacobian.setOldSolution(this->uOldTimeStep);
		NewtonMethod<G, ThisType> newtonMethod(this->grid, *this);
		newtonMethod.execute();
		*(this->uOldTimeStep) = *(this->u);
		
		return;
	}

  };

}
#endif
