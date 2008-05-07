#ifndef LEVELHETPROBLEM_HH
#define LEVELHETPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"
#include "dumux/material/randompermeability.hh"
#include <dune/istl/bvector.hh>

namespace Dune
{
	//! \ingroup diffusionProblems
	//! example class for diffusion problems
	template<class G, class RT>
	class LevelHetProblem : public DiffusionProblem<G,RT>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	
	public:
	  LevelHetProblem(G& g, const int level, const char* name = "permeab.dat", const bool create = true, 
			  				TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false)
	    : grid(g), DiffusionProblem<G,RT>(law, cap), permeability(g, level, name, create)
	  { }
	
	  const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi) 
	  {

		  return permeability.K(e);
	  }
	  
	  RT sat (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi)
	  {
		  return (*saturation)[grid.levelIndexSet(e.level()).index(e)];
	  }
	
	  RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi)
	  {
		return 0;
	  }
	
	  typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e, 
						   const Dune::FieldVector<DT,n>& xi) const
	  {
	    if ( x[0] > 300-1e-6 || x[0] < 1e-6) 
	    return Dune::BoundaryConditions::dirichlet;
	    // all other boundaries
	    return Dune::BoundaryConditions::neumann;
	  }
	
	  RT g (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		  return (x[0] < 1e-6) ? 2e6 : 1e6;
	  }
		  
	  RT J (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		if (x[0]<1e-6) return 1e-4;  
		return 0;
	  }
		  
		LevelRandomPermeability<G> permeability;
	private:
		G& grid;
	public:
		Dune::BlockVector<Dune::FieldVector<RT,1> >* saturation;
	};
}

#endif
