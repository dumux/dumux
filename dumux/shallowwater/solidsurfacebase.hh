// $Id: property_baseclasses.hh 703 2008-10-16 13:27:36Z jfritz $

#ifndef SOLIDSURFACEBASE_HH
#define SOLIDSURFACEBASE_HH

#include <dune/common/fvector.hh>
#include <vector>


namespace Dune
{

template<class G, class DT>
class SolidSurfaceBase
{
public:
	enum {dim=G::dimension};
	
	typedef typename G::Traits::template Codim<0>::Entity Entity;

	virtual DT frictionRelationType (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
	{
		DUNE_THROW(NotImplemented, "friction term not implemented!");
	}

	virtual DT friction (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
	{
		DUNE_THROW(NotImplemented, "friction term not implemented!");
	}

	virtual DT evalBottomElevation(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) = 0;
	
	virtual FieldVector<DT, dim> calcBottomSlopes(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) = 0;

	virtual ~SolidSurfaceBase()
	{}
};

} // end namespace
#endif /*PROPERTY_BASECLASSES*/
