// $Id$

#ifndef UNIFORMPROBLEM_HH
#define UNIFORMPROBLEM_HH

#include "dumux/diffusion/diffusionproblem.hh"
#include <dune/istl/bvector.hh>

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
	template<class G, class RT, class VC>
	class UniformProblem : public DiffusionProblem<G,RT,VC>
	{
	  typedef typename G::ctype DT;
	  enum {n=G::dimension};
	  typedef typename G::Traits::template Codim<0>::Entity Entity;

	public:
	  UniformProblem(VC& variableobj, TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false)
	    : DiffusionProblem<G,RT,VC>(variableobj,law, cap)
	    {
		permloc = 0;
	    for (int k = 0; k < n; k++)
	      permloc[k][k] = 1e-10;

	    }

	  UniformProblem()
	    : DiffusionProblem<G,RT,VC>()
	  {
	    permloc = 0;
	    for (int k = 0; k < n; k++)
	      permloc[k][k] = 1e-10;
	  }

	  Dune::FieldMatrix<DT,n,n>& K (const Dune::FieldVector<DT,n>& x, const Entity& e,
					  const Dune::FieldVector<DT,n>& xi)
	  {
		  return permloc;
	  }


	  RT source   (const Dune::FieldVector<DT,n>& x, const Entity& e,
					  const Dune::FieldVector<DT,n>& xi)
	  {
//		Dune::FieldVector<DT,n> m(150);
//		if ((x-m).two_norm()<1) return 1e-6;
		return 0;
	  }

	  typename Dune::BoundaryConditions::Flags bctype (const Dune::FieldVector<DT,n>& x, const Entity& e,
						   const Dune::FieldVector<DT,n>& xi) const
	  {
	    if (x[0] > 1-1E-6 || x[0] < 1e-6)
	      return Dune::BoundaryConditions::dirichlet;
	    // all other boundaries
	    return Dune::BoundaryConditions::neumann;
	  }

	  RT dirichletPress (const Dune::FieldVector<DT,n>& x, const Entity& e,
					const Dune::FieldVector<DT,n>& xi) const
	  {
		  return (x[0] < 1e-6) ? 2e5 : 1e5;
	  }

		RT dirichletSat (const FieldVector<DT,n>& x, const Entity& e,
			   const FieldVector<DT,n>& xi) const
		{
			if (x[0] < 1e-6)
				return 0.8;
			else
				return 0.2;
		}

	  RT neumannPress (const Dune::FieldVector<DT,n>& x, const Entity& e,
					const Dune::FieldVector<DT,n>& xi) const
	  {
		return 0;
	  }


	private:
		Dune::FieldMatrix<DT,n,n> permloc;
	};
}

#endif
