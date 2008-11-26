#ifndef TESTPROBLEM_UPSS_HH
#define TESTPROBLEM_UPSS_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"


namespace Dune {
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class G, class RT, class VC> class UpsSProblem :
	public FractionalFlowProblem<G,RT,VC> {

	typedef typename G::ctype DT;
	enum {n=G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
	UpsSProblem(VC& variableobj, Fluid& wp, Fluid& nwp, Matrix2p<G, RT>& s, TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G,RT>), const FieldVector<DT,n> Left = 0,
			const FieldVector<DT,n> Right = 0, const bool cap = false) :
		FractionalFlowProblem<G, RT, VC>(variableobj, wp, nwp, s, law, cap), Left_(Left[0]),
				Right_(Right[0]), eps_(1e-8)
				{}

	virtual RT qPress  (const FieldVector<DT,n>& x, const Entity& e,
							const FieldVector<DT,n>& xi)
	{
		return 0;
	}

	typename BoundaryConditions::Flags bctypePress(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) const {
		if ((x[0] < eps_ ) || x[0] > (Right_ - eps_))
			return BoundaryConditions::dirichlet;
		// all other boundaries
		return BoundaryConditions::neumann;
	}

    BoundaryConditions::Flags bctypeSat (const FieldVector<DT,n>& x, const Entity& e,
				      const FieldVector<DT,n>& xi) const
    {
      if (x[0] < eps_)
	return Dune::BoundaryConditions::dirichlet;
      else
	return Dune::BoundaryConditions::neumann;
    }

	RT gPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const {
		if (x[0] < eps_)
			return 1;
		// all other boundaries
		return 0;
	}

	RT gSat(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const {
		if (x[0] < eps_)
			return 1;
		// all other boundaries
		return 0;
	}

	RT JPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const {
//		if (x[0] > Right_ - eps_)
//			return 3e-7;
		return 0;
	}

	RT JSat(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi, RT& factor) const
			{
		if (x[0] > (Right_ - eps_))
			return factor;
		// all other boundaries
		return 0;
	}

    RT initSat (const FieldVector<DT,n>& x, const Entity& e,
	   const FieldVector<DT,n>& xi) const
    {
//      if (x[0] < eps_)
//	return 1;
//      else
	return 0;
    }

private:
	DT Left_;
	DT Right_;

	RT eps_;
};
}

#endif
