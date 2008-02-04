#ifndef TESTPROBLEM_HH
#define TESTPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
	template<class G, class RT>
	class TestProblem : public DiffusionProblem<G,RT>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	
	public:
	  TestProblem(G& g, TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false)
	    : DiffusionProblem<G,RT>(law, cap)
	  { }
	
	  TestProblem()
	    : DiffusionProblem<G,RT>()
	  { }
	
	  const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi) 
	  {
	    if (x[0] < 300) 
	      permloc[0][0] = permloc[1][1] = 1e-10;
	    else 
	      permloc[0][0] = permloc[1][1] = 2e-10;

	    permloc[0][1] = permloc[1][0] = 0;
	    
	    return permloc;
	  }
	
	  RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi) const
	  {
		return 0;
	  }
	
	  typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e, 
						   const Dune::FieldVector<DT,n>& xi) const
	  {
	    if (x[0] > 600-1E-6 || x[0] < 1e-6) 
	      return Dune::BoundaryConditions::dirichlet;
	    // all other boundaries
	    return Dune::BoundaryConditions::neumann;
	  }
	
	  RT g (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		  return (x[0] < 1e-6) ? 1e6 : 0;
	  }
		  
		
	  RT J (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		return 0;
	  }
		  
	private:
		Dune::FieldMatrix<DT,n,n> permloc;
	};
}

#endif
