#ifndef FVCA5TEST6PROBLEM_HH
#define FVCA5TEST6PROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
	template<class G, class RT, class VC>
	class FVCA5Test6Problem : public DiffusionProblem<G,RT,VC>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	
	public:
	  FVCA5Test6Problem(VC& variables, double delta = 0.2)
	    : DiffusionProblem<G,RT,VC>(variables)
	  { 
		  delta_ = delta;
		  theta_ = atan(delta);
	  }
	
	  FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
					  const FieldVector<DT,n>& xi) 
	  {
		  double phi1 = x[1] - delta_*(x[0] - 0.5) - 0.475; 
		  double phi2 = phi1 - 0.05; 
		  double alpha;
		  double beta;
		  if (phi1 < 0.0  || phi2 > 0.0) { 
		     alpha = 1.0; 
		     beta = 0.1;
		  }
		  else {
		     alpha = 100.0;
		     beta = 10.0;
		  }
		  double cost = 1.0/sqrt(1.0 + delta_*delta_);
		  double sint = delta_*cost;
		  permloc_[0][0] = alpha*cost*cost + beta*sint*sint;
		  permloc_[0][1] = permloc_[1][0] = cost*sint*(alpha - beta);
		  permloc_[1][1] = alpha*sint*sint + beta*cost*cost;

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
	      return BoundaryConditions::dirichlet;
	  }
	
	  RT g (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	  {
		  return (exact(x));
	  }
		  
	  RT J (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	  {
		return 0;
	  }

	  RT exact (const FieldVector<DT,n>& x) const
	  {
		  return (- x[0] - delta_*x[1]);
	  }

	  FieldVector<RT,n> exactGrad (const FieldVector<DT,n>& x) const
	  {	
		  FieldVector<RT,n> grad(0);
		  
		  grad[0] = -1.0;
		  grad[1] = -delta_;
		  
		  return grad;
	  }
	
	private:
		FieldMatrix<DT,n,n> permloc_;
		double delta_;
		double theta_;
	};
}

#endif
