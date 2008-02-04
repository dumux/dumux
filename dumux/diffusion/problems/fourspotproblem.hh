#ifndef FOURSPOTPROBLEM_HH
#define FOURSPOTPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"
#include "dumux/material/randompermeability.hh"

namespace Dune
{
	//! \ingroup diffusionProblems
	//! example class for diffusion problems
	template<class G, class RT>
	class FourSpotProblem : public DiffusionProblem<G,RT>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  
	  Dune::FieldMatrix<DT,n,n> K_;
	
	public:
	  FourSpotProblem(G& g, const int level, const char* name = "permeab.dat", const bool create = true, 
			  				TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false)
	    : DiffusionProblem<G,RT>(law, cap), permeability(g, level, name, create)
	  { }
	
	  const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi) 
	  {
//		  K_[0][0] = K_[1][1] = 1e-10;
//		  Dune::FieldVector<DT,n> center(155);
//		  double radius = 50;
//		  if ((x-center).two_norm()<radius) K_[0][0] = K_[1][1] = 1e-6;
//		  if (x[0]>150 && x[0]<200 && x[1]<145 || x[1]>150 && x[1]<200 && x[0]<145) K_[0][0] = K_[1][1] = 1e-14;
//		  return K_;
		  return permeability.K(e);
	  }
	
	  RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi)
	  {
		return 0;
	  }
	
	  typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e, 
						   const Dune::FieldVector<DT,n>& xi) const
	  {
	    if ( x[1]>285 && x[0]<15 || x[1]<15 && x[0]>285) return Dune::BoundaryConditions::dirichlet;
	    return Dune::BoundaryConditions::neumann;
	  }
	
	  RT g (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
//		  if (x[0]<15 && x[1]<15 || x[0]>285 && x[1]>285) return 1e7;
		  return  1e6;
	  }
		  
	  RT J (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		if (x[0]<15 && x[1]<15 || x[0]>285 && x[1]>285) return 1e-3;
		return 0;
	  }
		  
		LevelRandomPermeability<G> permeability;
	private:
	};
}

#endif
