#ifndef FVCA5TEST9PROBLEM_HH
#define FVCA5TEST9PROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
	template<class G, class RT>
	class FVCA5Test9Problem : public DiffusionProblem<G,RT>
	{
		template<int dim>
	    struct ElementLayout
	    {
	      bool contains (Dune::GeometryType gt)
	      {
	    	  return gt.dim() == dim;
	      }
	    }; 

	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	  typedef typename G::Traits::LevelIndexSet IS;
	  typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
	
	public:
	  FVCA5Test9Problem(double delta = 1.0e-3, double theta = 1.178097245096172418)
	    : DiffusionProblem<G,RT>()
	  { 
		  double cost = cos(theta); 
		  double sint = sqrt(1.0 - cost*cost); 

		  permloc_[0][0] = cost*cost + delta*sint*sint; 
		  permloc_[1][1] = sint*sint + delta*cost*cost;
		  permloc_[0][1] = permloc_[1][0] = cost*sint*(1.0 - delta);
	  }
	
	  const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
					  const FieldVector<DT,n>& xi) 
	  {
		  return permloc_;
	  }
	
	  RT q   (const FieldVector<DT,n>& x, const Entity& e, 
					  const FieldVector<DT,n>& xi)
	  {
		  return (0.0);  
	  }
	
	  typename BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e, 
						   const FieldVector<DT,n>& xi) const
	  {
		  if (fabs(x[0]) < 1e-8 || fabs(x[0] - 1.0) < 1e-8 
				  || fabs(x[1]) < 1e-8 || fabs(x[1] - 1.0) < 1e-8)
			  return BoundaryConditions::neumann;
		  
		  return BoundaryConditions::dirichlet;
	  }
	
	  RT g (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	  {
		  if (x[0] < 0.5)
			  return (0.0);
		  else 
			  return (1.0);
	  }
		  
		
	  RT J (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	  {
		return 0;
	  }

	  RT exact (const FieldVector<DT,n>& x) const
	  {
		  return 0;
	  }

	  FieldVector<RT,n> exactGrad (const FieldVector<DT,n>& x) const
	  {	
		  FieldVector<RT,n> grad(0);
		  return grad;
	  }
	
	private:
		FieldMatrix<DT,n,n> permloc_;
	};
}

#endif
