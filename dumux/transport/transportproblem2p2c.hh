#ifndef TRANSPORTPROBLEM2P2C_HH
#define TRANSPORTPROBLEM2P2C_HH

#include <iostream>
#include <iomanip>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/disc/operators/boundaryconditions.hh>
#include <dumux/material/twophaserelations.hh>
#include <dumux/material/linearlaw.hh>

#include <dumux/transport/transportproblem.hh>
#include <dumux/material/properties.hh>
#include <dumux/operators/boundaryconditions2p2c.hh>

namespace Dune
{

  template<class G, class RT>
  class TransportProblem2p2c
  {

	typedef typename G::ctype DT;
	enum {n=G::dimension, m=1, blocksize=2*G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	
  public:
	virtual BoundaryConditions2p2c::Flags cbctype (const FieldVector<DT,n>& x, const Entity& e, 
					   const FieldVector<DT,n>& xi) const = 0;
	
	virtual BoundaryConditions2p2c::Flags ictype (const FieldVector<DT,n>& x, const Entity& e, 
					   const FieldVector<DT,n>& xi) const = 0;
	
	virtual BoundaryConditions::Flags pbctype (const FieldVector<DT,n>& x, const Entity& e, 
					   const FieldVector<DT,n>& xi) const = 0;

	virtual const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
						const FieldVector<DT,n>& xi) = 0;
	
	virtual RT gZ (const FieldVector<DT,n>& x, const Entity& e, 
		   const FieldVector<DT,n>& xi) const = 0; 
	
	virtual RT gS (const FieldVector<DT,n>& x, const Entity& e, 
			   const FieldVector<DT,n>& xi) const = 0; 
	
	virtual RT gPress (const FieldVector<DT,n>& x, const Entity& e, 
		   const FieldVector<DT,n>& xi) const = 0; 
	
	virtual FieldVector<RT,2> J (const FieldVector<DT,n>& x, const Entity& e, 
		   const FieldVector<DT,n>& xi) const = 0;
	
	virtual FieldVector<RT,2> q (const FieldVector<DT,n>& x, const Entity& e, 
			   const FieldVector<DT,n>& xi) const = 0;
	  
	virtual RT S0 (const FieldVector<DT,n>& x, const Entity& e, 
			const FieldVector<DT,n>& xi) const;
	
	virtual RT Z1_0 (const FieldVector<DT,n>& x, const Entity& e, 
				const FieldVector<DT,n>& xi) const
        {
	  return 0;
	}
		
	virtual RT porosity ()
	{
		return 1.0;
	}
	
	const FieldVector<DT,n>& gravity()
	{
		return gravity_;
	}

	TransportProblem2p2c(Henry& h, TwoPhaseRelations& law = *(new LinearLaw), 
								 const bool cap = false) 
	: capillary(cap), materialLaw(law), henry(h)
	{	
	}
	
	virtual ~TransportProblem2p2c () {}
	
	FieldVector<DT,n> gravity_;
	Henry& henry;
	const bool capillary;
	TwoPhaseRelations& materialLaw;
  };

}
#endif
