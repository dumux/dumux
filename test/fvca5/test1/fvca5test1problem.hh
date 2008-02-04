#ifndef FVCA5TEST1PROBLEM_HH
#define FVCA5TEST1PROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
	template<class G, class RT>
	class FVCA5Test1Problem : public DiffusionProblem<G,RT>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	
	public:
	  FVCA5Test1Problem()
	    : DiffusionProblem<G,RT>()
	  { }
	
	  const Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi) 
	  {
		  permloc[0][0] = permloc[1][1] = 1.5;
		  permloc[0][1] = permloc[1][0] = 0.5;
		  
		  return permloc;
	  }
	
	  RT q   (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					  const Dune::FieldVector<DT,n>& xi)
	  {
		  double uxx = -2.0*x[1]*(1.0 - x[1])*16.0;
		  double uxy = (-2.0*x[0] + 1.0)*(-2.0*x[1] + 1.0)*16.0; 
		  double uyy = -2.0*x[0]*(1.0 - x[0])*16.0; 

		  double kxx = 1.5; 
		  double kxy = 0.5; 
		  double kyy = 1.5; 

		  return (-(kxx*uxx + 2.0*kxy*uxy + kyy*uyy));  
	  }
	
	  typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e, 
						   const Dune::FieldVector<DT,n>& xi) const
	  {
	      return Dune::BoundaryConditions::dirichlet;
	  }
	
	  RT g (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		  return (exact(x));
	  }
		  
		
	  RT J (const Dune::FieldVector<DT,n>& x, const Entity& e, 
					const Dune::FieldVector<DT,n>& xi) const
	  {
		return 0;
	  }
		  
	  RT exact (const Dune::FieldVector<DT,n>& x) const
	  {
		  return (16.0*x[0]*(1.0 - x[0])*x[1]*(1.0 - x[1]));
	  }

	private:
		Dune::FieldMatrix<DT,n,n> permloc;
	};
}

#endif
