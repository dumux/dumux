// $Id$
#ifndef DISPERSIONSUBPROBLEM_HH
#define DISPERSIONSUBPROBLEM_HH

#include "dumux/upscaledsaturation/preprocess/fractionalflowproblemsubprobs.hh"
#include "dumux/transport/fv/numericalflux.hh"

namespace Dune
{
//! \ingroup diffusionProblems
//! example class for diffusion problems
template<class G, class RT, class VC> class DiffSubProblemX: public FractionalFlowProblemSubProbs<
		G, RT, VC>
{

typedef	typename G::ctype DT;
	enum
	{	n=G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
	DiffSubProblemX(VC& variableobj, Fluid& wp, Fluid& nwp, Matrix2p<typename G::HostGridType, RT>& s, TwoPhaseRelations<typename G::HostGridType, RT>& law = *(new TwoPhaseRelations<G,RT>), FieldVector<DT,n> Left = 0,
			FieldVector<DT,n> Right = 0) :
	FractionalFlowProblemSubProbs<G, RT, VC>(variableobj, wp, nwp, s, law), LowerLeft_(Left),
	UpperRight_(Right), eps_(1e-6)
	{}

	virtual RT qPress (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi)
	{
		return 0;
	}

	typename BoundaryConditions::Flags bctypePress(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) const
	{
		if (x[0] < LowerLeft_[0] +  eps_ || x[0] > (UpperRight_[0] - eps_))
		return BoundaryConditions::dirichlet;
		// all other boundaries
		return BoundaryConditions::neumann;
	}

	BoundaryConditions::Flags bctypeSat (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[0] < LowerLeft_[0] + eps_)
		return Dune::BoundaryConditions::dirichlet;
		// all other boundaries
		return Dune::BoundaryConditions::neumann;
	}

	RT gPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[0] <  LowerLeft_[0] + eps_)
		return 1;
		// all other boundaries
		return 0;
	}

	RT gSat(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[0] <  LowerLeft_[0] + eps_)
		return 1;
		// all other boundaries
		return 0;
	}

	RT JPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 0;
	}

	RT JSat(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi, RT& factor) const
	{
		if (x[0] > (UpperRight_[0] - eps_))
			return factor;
		// all other boundaries
		return 0;
	}

	RT initSat (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 0;
	}

private:
	FieldVector<DT,n> LowerLeft_;
	FieldVector<DT,n> UpperRight_;

	RT eps_;
};

template<class G, class RT, class VC> class DiffSubProblemY: public FractionalFlowProblemSubProbs<
		G , RT, VC>
{

typedef	typename G::ctype DT;
	enum
	{	n=G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
	DiffSubProblemY(VC& variableobj, Fluid& wp, Fluid& nwp, Matrix2p<typename G::HostGridType, RT>& s, TwoPhaseRelations<typename G::HostGridType, RT>& law = *(new TwoPhaseRelations<G,RT>), FieldVector<DT,n> Left = 0,
			FieldVector<DT,n> Right = 0) :
	FractionalFlowProblemSubProbs<G, RT, VC>(variableobj, wp, nwp, s, law), LowerLeft_(Left),
	UpperRight_(Right), eps_(1e-8)
	{}

	virtual RT qPress (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi)
	{
		return 0;
	}

	typename BoundaryConditions::Flags bctypePress(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) const
	{
		if ((x[1] <  LowerLeft_[1] + eps_ || x[1] > (UpperRight_[1] - eps_)))
		return BoundaryConditions::dirichlet;
		// all other boundaries
		return BoundaryConditions::neumann;
	}

	BoundaryConditions::Flags bctypeSat (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[1] <  LowerLeft_[1] + eps_)
		return Dune::BoundaryConditions::dirichlet;
		else
		return Dune::BoundaryConditions::neumann;
	}

	RT gPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[1] < LowerLeft_[1] + eps_)
		return 1;
		// all other boundaries
		return 0;
	}

	RT gSat(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[1] < LowerLeft_[1] + eps_)
		return 1;
		// all other boundaries
		return 0;
	}

	RT JPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 0;
	}

	RT JSat(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi, RT& factor) const
	{
		if (x[1] > (UpperRight_[1] - eps_))
			return factor;
		// all other boundaries
		return 0;
	}

	RT initSat (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 0;
	}

private:
	FieldVector<DT,n> LowerLeft_;
	FieldVector<DT,n> UpperRight_;

	RT eps_;
};
template<class G, class RT, class VC> class DiffSubProblemXInt2: public FractionalFlowProblemSubProbs<
		G, RT, VC>
{

typedef	typename G::ctype DT;
	enum
	{	n=G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
	DiffSubProblemXInt2(VC& variableobj, Fluid& wp, Fluid& nwp, Matrix2p<typename G::HostGridType, RT>& s, TwoPhaseRelations<typename G::HostGridType, RT>& law = *(new TwoPhaseRelations<G,RT>), FieldVector<DT,n> Left = 0,
			FieldVector<DT,n> Right = 0) :
	FractionalFlowProblemSubProbs<G, RT, VC>(variableobj, wp, nwp, s, law), LowerLeft_(Left),
	UpperRight_(Right), eps_(1e-6)
	{}

	virtual RT qPress (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi)
	{
		return 0;
	}

	typename BoundaryConditions::Flags bctypePress(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) const
	{
		if (x[0] < LowerLeft_[0] +  eps_ || x[0] > (UpperRight_[0] - eps_))
		return BoundaryConditions::dirichlet;
		// all other boundaries
		return BoundaryConditions::neumann;
	}

	BoundaryConditions::Flags bctypeSat (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[0] > UpperRight_[0] - eps_)
		return Dune::BoundaryConditions::dirichlet;
		// all other boundaries
		return Dune::BoundaryConditions::neumann;
	}

	RT gPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[0] <  LowerLeft_[0] + eps_)
		return 0;
		// all other boundaries
		return 1;
	}

	RT gSat(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[0] <  LowerLeft_[0] + eps_)
		return 0;
		// all other boundaries
		return 1;
	}

	RT JPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 0;
	}

	RT JSat(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi, RT& factor) const
	{
		if (x[0] < (LowerLeft_[0] + eps_))
			return factor;
		// all other boundaries
		return 0;
	}

	RT initSat (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 0;
	}

private:
	FieldVector<DT,n> LowerLeft_;
	FieldVector<DT,n> UpperRight_;

	RT eps_;
};

template<class G, class RT, class VC> class DiffSubProblemYInt2: public FractionalFlowProblemSubProbs<
		G , RT, VC>
{

typedef	typename G::ctype DT;
	enum
	{	n=G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
	DiffSubProblemYInt2(VC& variableobj, Fluid& wp, Fluid& nwp, Matrix2p<typename G::HostGridType, RT>& s, TwoPhaseRelations<typename G::HostGridType, RT>& law = *(new TwoPhaseRelations<G,RT>), FieldVector<DT,n> Left = 0,
			FieldVector<DT,n> Right = 0) :
	FractionalFlowProblemSubProbs<G, RT, VC>(variableobj, wp, nwp, s, law), LowerLeft_(Left),
	UpperRight_(Right), eps_(1e-8)
	{}

	virtual RT qPress (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi)
	{
		return 0;
	}

	typename BoundaryConditions::Flags bctypePress(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) const
	{
		if ((x[1] <  LowerLeft_[1] + eps_ || x[1] > (UpperRight_[1] - eps_)))
		return BoundaryConditions::dirichlet;
		// all other boundaries
		return BoundaryConditions::neumann;
	}

	BoundaryConditions::Flags bctypeSat (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[1] >  UpperRight_[1] - eps_)
		return Dune::BoundaryConditions::dirichlet;
		else
		return Dune::BoundaryConditions::neumann;
	}

	RT gPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[1] < LowerLeft_[1] + eps_)
		return 0;
		// all other boundaries
		return 1;
	}

	RT gSat(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		if (x[1] < LowerLeft_[1] + eps_)
		return 0;
		// all other boundaries
		return 1;
	}

	RT JPress(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 0;
	}

	RT JSat(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi, RT& factor) const
	{
		if (x[1] < (LowerLeft_[1] + eps_))
			return factor;
		// all other boundaries
		return 0;
	}

	RT initSat (const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const
	{
		return 0;
	}

private:
	FieldVector<DT,n> LowerLeft_;
	FieldVector<DT,n> UpperRight_;

	RT eps_;
};
}
#endif
