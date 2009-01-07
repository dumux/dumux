#ifndef FVCA5TEST4PROBLEM_HH
#define FVCA5TEST4PROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
	template<class G, class RT, class VC>
	class FVCA5Test4Problem : public DiffusionProblem<G,RT,VC>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;

	public:
	  FVCA5Test4Problem(VC& variables)
	    : DiffusionProblem<G,RT,VC>(variables)
	  { }

	  FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e,
					  const FieldVector<DT,n>& xi)
	  {
		  if (x[0] < 0.5) {
		     int inty = (int)(10.0*(x[1] + 0.15));
		     if (inty - 2*(int)(inty/2) == 0) {
		        permloc_[0][0] = 100.0;
		        permloc_[0][1] = permloc_[1][0] = 0.0;
		        permloc_[1][1] = 10.0;
		     }
		     else {
		        permloc_[0][0] = 0.01;
		        permloc_[0][1] = permloc_[1][0] = 0.0;
		        permloc_[1][1] = 0.001;
		     }
		  }
		  else {
		     int inty = (int)(10.0*x[1]);
		     if (inty - 2*(int)(inty/2) == 0) {
		        permloc_[0][0] = 100.0;
		        permloc_[0][1] = permloc_[1][0] = 0.0;
		        permloc_[1][1] = 10.0;
		     }
		     else {
		        permloc_[0][0] = 0.01;
		        permloc_[0][1] = permloc_[1][0] = 0.0;
		        permloc_[1][1] = 0.001;
		     }
		  }

		  return permloc_;
	  }

	  RT source   (const FieldVector<DT,n>& x, const Entity& e,
					  const FieldVector<DT,n>& xi)
	  {
		  return (0.0);
	  }

	  typename BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e,
						   const FieldVector<DT,n>& xi) const
	  {
	      return BoundaryConditions::dirichlet;
	  }

	  RT dirichletPress (const FieldVector<DT,n>& x, const Entity& e,
					const FieldVector<DT,n>& xi) const
	  {
		  return (1.0 - x[0]);
	  }

	  RT neumannPress (const FieldVector<DT,n>& x, const Entity& e,
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
