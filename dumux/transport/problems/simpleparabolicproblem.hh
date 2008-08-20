#ifndef DUNE_SIMPLEPARABOLICPROBLEM_HH
#define DUNE_SIMPLEPARABOLICPROBLEM_HH

#include "dumux/transport/transportproblem.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief example class for a transport problem
  template<class G, class RT, class VC>
  class SimpleParabolicProblem : public TransportProblem<G, RT, VC > {
	  typedef typename G::ctype DT;
	  enum {n=G::dimension, m=1};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef typename Dune::BlockVector< Dune::FieldVector<Dune::FieldVector<double, n>, 2*n> > VelType;
  private:
	  DT left;
	  DT right;
	  FieldVector<DT,n> vLoc;

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
	  
	RT S0 (const FieldVector<DT,n>& x, const Entity& e, 
			const FieldVector<DT,n>& xi) const 
	{
		return 0;
	}
	  
	const FieldVector<DT,n>& vTotal (const Entity& e, const int numberInSelf)
	{
		vLoc[0] = 0;
		vLoc[1] = 0;
		
		return vLoc;
	}

	SimpleParabolicProblem(VC& variables, const G& g, TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false) 
	: TransportProblem<G, RT, VC>(variables, law, cap), left((g.lowerLeft())[0]), right((g.upperRight())[0])
	{	}

	SimpleParabolicProblem(VC& variables, TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false) 
	: TransportProblem<G, RT, VC>(variables, law, cap), left(0), right(1)
	{	}
  };

}
#endif
