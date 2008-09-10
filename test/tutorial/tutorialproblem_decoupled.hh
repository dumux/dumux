#ifndef TUTORIALPROBLEM_DECOUPLED_HH
#define TUTORIALPROBLEM_DECOUPLED_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class G, class RT, class VC> class TutorialProblemDecoupled: public FractionalFlowProblem<
		G, RT, VC>
{

typedef	typename G::ctype DT;
	enum
	{	n=G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
	TutorialProblemDecoupled(VC& variableobj, MediumNonIsothermal& wp, MediumNonIsothermal& nwp, Matrix2p<G, RT>& s, TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G,RT>), const FieldVector<DT,n> Left = 0,
			const FieldVector<DT,n> Right = 0) :
	FractionalFlowProblem<G, RT, VC>(variableobj, wp, nwp, s, law), Left_(Left[0]),
	Right_(Right[0]), eps_(1e-8)
	{}

	virtual RT qPress (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi)
	{
		return 0;
	}

	typename BoundaryConditions::Flags bctypePress(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) const
	{
		if ((x[0] < eps_))
		return BoundaryConditions::dirichlet;
		// all other boundaries
		return BoundaryConditions::neumann;
	}

	BoundaryConditions::Flags bctypeSat (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[0]> (Right_ - eps_) || x[0] < eps_)
		return Dune::BoundaryConditions::dirichlet;
		else
		return Dune::BoundaryConditions::neumann;
	}

	RT gPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 2e5;
	}

	RT gSat(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[0] < eps_)
		return 1;
		// all other boundaries
		return 0;
	}

	RT JPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[0]> Right_ - eps_)
		return 3e-7;
		// all other boundaries
		return 0;
	}

	RT S0 (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 0.2;
	}

private:
	DT Left_;
	DT Right_;

	RT eps_;
};
}

#endif
