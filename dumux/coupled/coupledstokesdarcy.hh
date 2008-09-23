// $Id$

#ifndef DUNE_COUPLEDSTOKESDARCY_HH
#define DUNE_COUPLEDSTOKESDARCY_HH

#include "coupledmodel.hh"

namespace Dune
{

template<class StokesModel, class DarcyModel>
class CoupledStokesDarcy : public CoupledModel<StokesModel, DarcyModel, CoupledStokesDarcy<StokesModel, DarcyModel> > {
public:
	// typedefs depending on base class 
	typedef CoupledModel<StokesModel, DarcyModel, CoupledStokesDarcy<StokesModel, DarcyModel> > BaseType;
	typedef typename BaseType::FirstGrid StokesGrid;
	typedef typename BaseType::SecondGrid DarcyGrid;

	template <class A12Type, class A21Type>
	void assembleCoupling(A12Type& A12, A21Type& A21)
	{
		// intialize rowsizes with 0
		for (typename A12Type::size_type i = 0; i < A12.N(); i++)
			A12.setrowsize(i, 0);
		for (typename A21Type::size_type i = 0; i < A21.N(); i++)
			A21.setrowsize(i, 0);
		A12.endrowsizes();
		A21.endrowsizes();
		A12.endindices();
		A21.endindices();
	}

	CoupledStokesDarcy(const StokesGrid& stokesGrid, StokesModel& stokesModel,
			const DarcyGrid& darcyGrid, DarcyModel& darcyModel,
			bool assembleGlobalSystem)
	: BaseType(stokesGrid, stokesModel, darcyGrid, darcyModel, assembleGlobalSystem)
	{}
};

}

#endif
