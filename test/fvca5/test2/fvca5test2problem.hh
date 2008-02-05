#ifndef FVCA5TEST2PROBLEM_HH
#define FVCA5TEST2PROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
	template<class G, class RT>
	class FVCA5Test2Problem : public DiffusionProblem<G,RT>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;
	
	public:
	  FVCA5Test2Problem(double delta = 1.0e-3)
	    : DiffusionProblem<G,RT>()
	  { 
		  delta_ = delta;
		  permloc_[1][1] = 1.0; 
		  permloc_[0][0] = delta;
		  permloc_[1][0] = permloc_[0][1] = 0.0;
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
//		  if (x[1] < 1e-8)
//			  return BoundaryConditions::dirichlet;
		  
		  return BoundaryConditions::neumann;
	  }
	
	  RT g (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	  {
		  return (exact(x));
	  }
		  
		
	  RT J (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const
	  {
		  FieldVector<RT,n> unitNormal(0);
		  
		  if (x[0] < 1e-8) 
			  unitNormal[0] = -1.0;
		  else if (fabs(x[0] - 1.0) < 1e-8) 
			  unitNormal[0] = 1.0;
		  else if (x[1]  < 1e-8) 
			  unitNormal[1] = -1.0;
		  else 
			  unitNormal[1] = 1.0;
		  
		  FieldVector<RT,n> Kgrad(0);
		  permloc_.umv(exactGrad(x), Kgrad);

		  return -(Kgrad*unitNormal);
	  }

	  RT exact (const FieldVector<DT,n>& x) const
	  {
		  double pi = 4.0*atan(1.0); 
		  double x1 = 2.0*pi;  
		  double x2 = x1/sqrt(delta_); 
		  
		  return (sin(x1*x[1])*exp(-x2*x[0])); 
	  }

	  FieldVector<RT,n> exactGrad (const FieldVector<DT,n>& x) const
	  {	
		  FieldVector<RT,n> grad(0);
		  double pi = 4.0*atan(1.0); 
		  double x1 = 2.0*pi;  
		  double x2 = x1/sqrt(delta_); 
		  double u = sin(x1*x[1])*exp(-x2*x[0]);
		  grad[1] = x1*cos(x1*x[1])*exp(-x2*x[0]); 
		  grad[0] = -x2*u;
		  
		  return grad;
	  }
	
	private:
		FieldMatrix<DT,n,n> permloc_;
		double delta_;
	};
}

#endif
