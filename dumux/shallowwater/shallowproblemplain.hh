// $Id: buckleyleverettproblem.hh 621 2008-09-22 08:22:56Z markus $

#ifndef DUNE_SHALLOWPROBLEMPLAIN_HH
#define DUNE_SHALLOWPROBLEMPLAIN_HH

#include "dumux/shallowwater/shallowproblembase.hh"
#include"dumux/shallowwater/shallowvariableclass.hh"

namespace Dune
{
//! \ingroup transportProblems
//! @brief example class for a transport problem in shallwo water
template<class G, class DT, class VC> class ShallowProblemPlain :
	public ShallowProblemBase<G, DT,VC>
{
	enum
	{	dim=G::dimension, m=1, blocksize=2*G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef Dune::FieldVector<DT, dim> VelType;
	typedef Dune::FieldVector<DT,dim+1> SystemType;
	typedef Dune::SolidSurfaceBase<G,DT> Surface;

private:
	FieldVector<DT,dim> lowerLeft_;
	FieldVector<DT,dim> upperRight_;
	DT eps;

public:
	BoundaryConditions::Flags bctype(const FieldVector<DT,dim>& x,
			const Entity& e, const FieldVector<DT,dim>& xi) const
	{
		
		if (x[0] < eps)
		{
			return Dune::BoundaryConditions::dirichlet;
		}
		else
		{
			return Dune::BoundaryConditions::neumann;
		}
	}

	DT dirichlet(const FieldVector<DT,dim>& x, const Entity& e,
			const FieldVector<DT,dim>& xi) const
	{
		return 0;
	}

	SystemType neumann(const FieldVector<DT,dim>& x, const Entity& e,
			const FieldVector<DT,dim>& xi) const
	{
		SystemType boundaryFlux;

		if (x[0]= upperRight_[0] - eps)
		{
			boundaryFlux[0]=0;
			boundaryFlux[1]=0;
			boundaryFlux[2]=0;
			
		}
		return boundaryFlux;
			
	}

	DT setInitWDepth(const FieldVector<DT,dim>& x, const Entity& e,
			const FieldVector<DT,dim>& xi) const
	{
		return 0.00;
	}

	VelType setInitVel(const FieldVector<DT,dim>& x, const Entity& e,
			const FieldVector<DT,dim>& xi) const
	{
		VelType initVel_(0);
		return initVel_;

	}

	DT setSource(const FieldVector<DT,dim>& x, const Entity& e,
			const FieldVector<DT,dim>& xi) const
	{
		return 125E-7; //der Umrechnungsfaktor von l/(s ha) auf m/s ist 1E-7, Standardwert für Stuttgart 125 l/(s ha)
	}

	ShallowProblemPlain(VC& variableobject, Surface& surfaceobject,
			FieldVector<DT,dim>& L, FieldVector<DT,dim>& H) :
		ShallowProblemBase<G, DT, VC>(variableobject, surfaceobject),
				lowerLeft_(L), upperRight_(H), eps(1e-8)
	{
	}

};

}
#endif
