#ifndef GRAVITYPROBLEM_HH
#define GRAVITYPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
	template<class G, class RT>
	class GravityProblem : public DiffusionProblem<G,RT>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	
	public:
	  GravityProblem(G& g, TwoPhaseRelations& law = *(new LinearLaw), FieldVector<DT,n> gravity = *(new FieldVector<DT,n>(0)))
	    : DiffusionProblem<G,RT>(law, false, gravity)
	  { }
	
	  GravityProblem()
	    : DiffusionProblem<G,RT>()
	  { }
	
	  const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
					  const FieldVector<DT,n>& xi) 
	  {
		  permloc[0][0] = permloc[1][1] = 1e-10;
		  permloc[0][1] = permloc[1][0] = 0;
		  
		  return permloc;
	  }
	
	  RT q   (const FieldVector<DT,n>& x, const Entity& e, 
						  const FieldVector<DT,n>& xi) 
	  {
			  return 0;
	  }
	
	  typename BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e, 
						   const FieldVector<DT,n>& xi) const
	  {
	    if (x[1] > 10-1E-6)// || x[1] < 1e-6)
	      return BoundaryConditions::dirichlet;
	    // all other boundaries
	    return BoundaryConditions::neumann;
	  }
	
	  RT g (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	  {
		  //if (x[1] > 10-1e-6)
			  return (2e5);
	  }
		  
		
	  RT J (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	  {
//		  if (x[1] < 1e-8) 
//			  return (((this->gravity_)[n-1]));
//		  else 
			  return 0;
	  }
		  
	private:
		FieldMatrix<DT,n,n> permloc;
	};
}

#endif
