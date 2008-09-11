// $Id$ 

#ifndef DUNE_BUCKLEYLEVERETTPROBLEM_HH
#define DUNE_BUCKLEYLEVERETTPROBLEM_HH

#include "dumux/transport/transportproblem.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief example class for a transport problem
  template<class G, class RT, class VC>
  class BuckleyLeverettProblem 
  : public TransportProblem<G, RT,VC> {
		  
	  typedef typename G::ctype DT;
	  enum {n=G::dimension, m=1, blocksize=2*G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef Dune::FieldVector<double, n> R1;

  private:
	  DT left;
	  DT right;

  public:
	BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e, 
					   const FieldVector<DT,n>& xi) const
	{
		if (x[0] > right-1E-8 || x[0] < left+1e-8) 
			return Dune::BoundaryConditions::dirichlet;
		else
			return Dune::BoundaryConditions::neumann;
	}

	RT g (const FieldVector<DT,n>& x, const Entity& e, 
		   const FieldVector<DT,n>& xi) const 
	{
		if (x[0] < left+1e-8) 
			return 0.8;
		else
			return 0.2;
	}
	  
	RT S0 (const FieldVector<DT,n>& x, const Entity& e, 
			const FieldVector<DT,n>& xi) const 
	{
		return 0.2;
	}
	  

	BuckleyLeverettProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw), 
								const int level = 0, const bool cap = false) 
	: TransportProblem<G, RT, VC>(variableobj,law, cap), left((variableobj.grid.lowerLeft())[0]), right((variableobj.grid.upperRight())[0])
	{}
  };

}
#endif
