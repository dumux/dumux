#ifndef BRINKMANTESTPROBLEM_HH
#define BRINKMANTESTPROBLEM_HH

#include "dumux/brinkman/brinkmanproblem.hh"

namespace Dune
{
//! \ingroup BrinkmanProblems
//! example class for Brinkman problems
	template<class G, class RT>
	class BrinkmanTestProblem : public BrinkmanProblem<G,RT>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	
	public:
	  BrinkmanTestProblem()
	    : BrinkmanProblem<G,RT>()
	  {
	    permloc = 0; 
	    for (int k = 0; k < n; k++)
	      permloc[k][k] = 1e0;
	  }
	
	  const Dune::FieldMatrix<DT,n,n>& Kinv (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi) 
	  {
		  return permloc;
	  }
	  
	  virtual RT muEff   (const FieldVector<DT,n>& x, const Entity& e, 
			  const FieldVector<DT,n>& xi) 
	  {
		  return 1.0;
	  }
	  
	  virtual RT mu   (const FieldVector<DT,n>& x, const Entity& e, 
						const FieldVector<DT,n>& xi) 
	  {
		  return 1.0;
	  }

	  RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi)
	  {
		return 0;
	  }
	
	  typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e, 
						   const Dune::FieldVector<DT,n>& xi) const
	  {
	    if (x[0] > 2-1E-6 || x[0] < 1e-6) 
	      return Dune::BoundaryConditions::dirichlet;

	    return Dune::BoundaryConditions::neumann;
	  }
	
	  RT g (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		  return (x[0] < 1e-6) ? 2e5 : 1e5;
	  }
		  
		
	  RT J (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		return 0;
	  }
		  
	  
	private:
		Dune::FieldMatrix<DT,n,n> permloc;
		G* grid;
	};
}

#endif
