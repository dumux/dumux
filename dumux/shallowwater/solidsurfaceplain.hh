// $Id: matrixproperties.hh 605 2008-09-18 16:43:57Z melanie $

#ifndef SOLIDSURFACEPLAIN_HH
#define SOLIDSURFACEPLAIN_HH

#include <dumux/shallowwater/solidsurfacebase.hh>

namespace Dune
{

template<class G, class DT> class SolidSurfacePlain :
	public SolidSurfaceBase<G,DT>
{
public:
	typedef typename G::Traits::template Codim<0>::Entity Entity;

	enum
	{	dim=G::dimension};

	DT evalBottomElevation(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
	{

		bottomElevation = 100 + crossSlope*x[0]+ longSlope*x[1];

		return bottomElevation;

	}
	FieldVector<DT, dim> calcBottomSlopes(const FieldVector<DT,dim>& x,
			const Entity& e, const FieldVector<DT,dim>& xi)
	{

		bottomSlope[0]=crossSlope;
		bottomSlope[1]=longSlope;

		return bottomSlope;
	}

	SolidSurfacePlain() :
		SolidSurfaceBase<G, DT>(), crossSlope(-0.02), longSlope(0.0),eps(1e-12)
	{

	}

private:
	//Definition der Variablen für Quer- und Längsneigung
	DT bottomElevation;
	DT crossSlope;
	DT longSlope;
	FieldVector<DT, dim> bottomSlope;
	DT eps;

};

} // end namespace
#endif 
