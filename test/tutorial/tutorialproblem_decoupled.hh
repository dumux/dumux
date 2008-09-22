#ifndef TUTORIALPROBLEM_DECOUPLED_HH
#define TUTORIALPROBLEM_DECOUPLED_HH

#include "dumux/fractionalflow/fractionalflowproblem.hh"

namespace Dune
{

template<class G, class RT, class VC> class TutorialProblemDecoupled /*@\label{tutorial-decoupled:tutorialproblem}@*/
	: public FractionalFlowProblem<G, RT, VC>
{

typedef	typename G::ctype DT;
	enum
	{	n=G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
	TutorialProblemDecoupled(VC& variableobj, Fluid& wp, Fluid& nwp, Matrix2p<G, RT>& s,
			TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G,RT>),
			const FieldVector<DT,n> Left = 0, const FieldVector<DT,n> Right = 0)
	: FractionalFlowProblem<G, RT, VC>(variableobj, wp, nwp, s, law),
	Left_(Left[0]), Right_(Right[0]), eps_(1e-8)
	{}

	// function returning source/sink terms for the pressure equation
	// depending on the position within the domain
	virtual RT qPress (const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-decoupled:qpress}@*/
			const FieldVector<DT,n>& xi)
	{
		return 0;
	}

	// function returning the boundary condition type for solution
	// of the pressure equation depending on the position within the domain
	typename BoundaryConditions::Flags bctypePress(const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-decoupled:bctypepress}@*/
			const FieldVector<DT,n>& xi) const
	{
		if (x[0] < eps_)
		{
			return BoundaryConditions::dirichlet;
		}
		// all other boundaries
		return BoundaryConditions::neumann;
	}

	// function returning the boundary condition type for solution
	// of the saturation equation depending on the position within the domain
	BoundaryConditions::Flags bctypeSat (const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-decoupled:bctypesat}@*/
			const FieldVector<DT,n>& xi) const
	{
		if (x[0]> (Right_ - eps_) || x[0] < eps_)
		{
			return Dune::BoundaryConditions::dirichlet;
		}
		// all other boundaries
		return Dune::BoundaryConditions::neumann;
	}

	// function returning the Dirichlet boundary condition for the solution
	// of the pressure equation depending on the position within the domain
	RT gPress(const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-decoupled:gpress}@*/
			const FieldVector<DT,n>& xi) const
	{
		return 2e5;
	}

	// function returning the Dirichlet boundary condition for the solution
	// of the saturation equation depending on the position within the domain
	RT gSat(const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-decoupled:gsat}@*/
			const FieldVector<DT,n>& xi) const
	{
		if (x[0] < eps_)
		{
			return 1;
		}
		// all other boundaries
		return 0;
	}

	// function returning the Neumann boundary condition for the solution
	// of the pressure equation depending on the position within the domain
	RT JPress(const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-decoupled:jpress}@*/
			const FieldVector<DT,n>& xi) const
	{
		if (x[0]> Right_ - eps_)
		{
			return 3e-7;
		}
		// all other boundaries
		return 0;
	}

	// function returning the initial saturation
	// depending on the position within the domain
	RT initSat (const FieldVector<DT,n>& x, const Entity& e, /*@\label{tutorial-decoupled:initsat}@*/
			const FieldVector<DT,n>& xi) const
	{
		return 0;
	}

private:
	DT Left_;
	DT Right_;

	RT eps_;
};
} // end namespace
#endif
