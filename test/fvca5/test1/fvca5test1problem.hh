#ifndef FVCA5TEST1PROBLEM_HH
#define FVCA5TEST1PROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
	template<class G, class RT, class VC>
	class FVCA5Test1Problem : public DiffusionProblem<G,RT,VC>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;

	public:
	  FVCA5Test1Problem(VC& variables)
	    : DiffusionProblem<G,RT,VC>(variables)
	  { }

	  virtual FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e,
					  const FieldVector<DT,n>& xi)
	  {
		  permloc[0][0] = permloc[1][1] = 1.5;
		  permloc[0][1] = permloc[1][0] = 0.5;

		  return permloc;
	  }

	  RT q   (const FieldVector<DT,n>& x, const Entity& e,
					  const FieldVector<DT,n>& xi)
	  {
		  // Test 1.1:
		  double uxx = -2.0*x[1]*(1.0 - x[1])*16.0;
		  double uxy = (-2.0*x[0] + 1.0)*(-2.0*x[1] + 1.0)*16.0;
		  double uyy = -2.0*x[0]*(1.0 - x[0])*16.0;

		  // Test 1.2:
//		  double x1 = 1.0 - x[0];
//		  double y1 = 1.0 - x[1];
//		  double uxx = -y1*y1*sin(x1*y1) + 6.0*x1*y1*y1;
//		  double uxy = -x1*y1*sin(x1*y1) + cos(x1*y1) + 6.0*y1*x1*x1;
//		  double uyy = -x1*x1*sin(x1*y1) + 2.0*x1*x1*x1;

		  double kxx = 1.5;
		  double kxy = 0.5;
		  double kyy = 1.5;

		  return (-(kxx*uxx + 2.0*kxy*uxy + kyy*uyy));
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

	  RT gSat (const FieldVector<DT,n>& x, const Entity& e,
					const FieldVector<DT,n>& xi) const
	  {
		  return (1);
	  }


	  RT J (const FieldVector<DT,n>& x, const Entity& e,
					const FieldVector<DT,n>& xi) const
	  {
		return 0;
	  }

	  RT exact (const FieldVector<DT,n>& x) const
	  {
		  // Test 1.1:
		  return (16.0*x[0]*(1.0 - x[0])*x[1]*(1.0 - x[1]));

		  // Test 1.2:
//		  double x1 = 1.0 - x[0];
//		  double y1 = 1.0 - x[1];
//		  return (sin(x1*y1) + x1*x1*x1*y1*y1);
	  }

	  FieldVector<RT,n> exactGrad (const FieldVector<DT,n>& x) const
	  {
		  FieldVector<RT,n> grad;

		  // Test 1.1:
		  grad[0] = 16.0*(1.0 - 2.0*x[0])*(x[1] - x[1]*x[1]);
		  grad[1] = 16.0*(1.0 - 2.0*x[1])*(x[0] - x[0]*x[0]);

		  // Test 1.2:
//		  double x1 = 1.0 - x[0];
//		  double y1 = 1.0 - x[1];
//		  grad[0] = -y1*cos(x1*y1) - 3.0*(x1*y1)*(x1*y1);
//		  grad[1] = -x1*cos(x1*y1) - 2.0*y1*x1*x1*x1;

		  return grad;
	  }

	private:
		FieldMatrix<DT,n,n> permloc;
	};
}

#endif
