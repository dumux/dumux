#ifndef UNIFORMPROBLEM_HH
#define UNIFORMPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
	template<class G, class RT>
	class UniformProblem : public DiffusionProblem<G,RT>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	
	public:
	  UniformProblem(G& g, TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false)
	    : DiffusionProblem<G,RT>(law, cap)
	  { }
	
	  UniformProblem()
	    : DiffusionProblem<G,RT>()
	  {
	    permloc = 0; 
	    for (int k = 0; k < n; k++)
	      permloc[k][k] = 1e-10;
	  }
	
	  const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi) 
	  {
		  return permloc;
	  }
	
	  RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi)
	  {
		return 0;
	  }
	
	  typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e, 
						   const Dune::FieldVector<DT,n>& xi) const
	  {
	    if (x[0] > 10-1E-6 || x[0] < 1e-6) 
	      return Dune::BoundaryConditions::dirichlet;
	    // all other boundaries
	    return Dune::BoundaryConditions::neumann;
	  }
	
	  RT g (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		  return (x[0] < 1e-6) ? 1.1e6 : 1e6;
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
