// $Id$

#ifndef DUNE_SIMPLEPROBLEM_HH
#define DUNE_SIMPLEPROBLEM_HH

#include "dumux/transport/transportproblem.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief example class for a transport problem
  template<class G, class RT, class VC>
  class SimpleProblem : public TransportProblem<G, RT, VC> {
	  typedef typename G::ctype DT;
	  enum {n=G::dimension, m=1};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;

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
			return 1;
		else
			return 0;
	}

	RT initSat (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 0;
	}

	SimpleProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw), FieldVector<DT,n>& Left = 0, FieldVector<DT,n>& Right = 1, const bool cap = false)
	: TransportProblem<G, RT, VC>(variableobj,law, cap), left(Left[0]), right(Right[0])
	{	}

	SimpleProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false)
	: TransportProblem<G, RT, VC>(variableobj,law, cap), left(0), right(1)
	{	}
  };

}
#endif
