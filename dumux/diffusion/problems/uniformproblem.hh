#ifndef UNIFORMPROBLEM_HH
#define UNIFORMPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"
#include <dune/istl/bvector.hh>

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
	  UniformProblem(G& g, const bool sflag = true, TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false)
	    : DiffusionProblem<G,RT>(law, cap), saturationflag(sflag), grid(&g) 
	  { 
		permloc = 0;
	    for (int k = 0; k < n; k++)
	      permloc[k][k] = 1e-10;
	    
	    }
	
	  UniformProblem()
	    : DiffusionProblem<G,RT>()
	  {
	    permloc = 0; 
	    for (int k = 0; k < n; k++)
	      permloc[k][k] = 1e-10;
	    
	    saturationflag = false;
	  }
	
	  const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi) 
	  {
		  return permloc;
	  }
	  
	  RT sat (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi)
	  {
		  if (saturationflag)
			  return (*saturation)[grid->levelIndexSet(e.level()).index(e)];
		  return 0;
	  }
	
	  RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi)
	  {
//		Dune::FieldVector<DT,n> m(150);
//		if ((x-m).two_norm()<1) return 1e-6;
		return 0;
	  }
	
	  typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e, 
						   const Dune::FieldVector<DT,n>& xi) const
	  {
	    if (x[0] > 300-1E-6 || x[0] < 1e-6) 
	      return Dune::BoundaryConditions::dirichlet;
	    // all other boundaries
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
		bool saturationflag;
		G* grid;
	public:
		Dune::BlockVector<Dune::FieldVector<RT,1> >* saturation;
	};
}

#endif
