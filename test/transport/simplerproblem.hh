#ifndef DUNE_SIMPLERPROBLEM_HH
#define DUNE_SIMPLERPROBLEM_HH

#include "dumux/transport/transportproblem.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief example class for a transport problem
  template<class G, class RT>
  class SimplerProblem : public TransportProblem<G, RT, Dune::BlockVector< Dune::FieldVector<Dune::FieldVector<double, G::dimension>, 2*G::dimension> > > {
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

	RT initSat (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 0;
	}

	const FieldVector<DT,n>& vTotal (const Entity& e, const int numberInSelf)
	{
		vLoc[0] = 1.0/6.0*1e-6;
		vLoc[1] = 0;

		return vLoc;
	}

	SimplerProblem(const G& g, TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false)
	: TransportProblem<G, RT, VelType>(law, cap), left(0), right(600)
	{	}

	SimplerProblem(TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false)
	: TransportProblem<G, RT, VelType>(law, cap), left(0), right(1)
	{	}
  };

}
#endif
