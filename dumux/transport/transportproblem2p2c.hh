#ifndef TRANSPORTPROBLEM2P2C_HH
#define TRANSPORTPROBLEM2P2C_HH

#include "dumux/transport/transportproblem.hh"
#include "dumux/material/properties.hh"
#include "dumux/operators/boundaryconditions2p2c.hh"

namespace Dune
{

  template<class G, class RT, class VelType>
  class TransportProblem2p2c : public TransportProblem<G, RT, VelType>
  {

	typedef typename G::ctype DT;
	enum {n=G::dimension, m=1, blocksize=2*G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	
  public:
	virtual BoundaryConditions2p2c::Flags cbctype (const FieldVector<DT,n>& x, const Entity& e, 
					   const FieldVector<DT,n>& xi) const = 0;
	
	virtual BoundaryConditions2p2c::Flags ictype (const FieldVector<DT,n>& x, const Entity& e, 
					   const FieldVector<DT,n>& xi) const = 0;

	virtual RT g (const FieldVector<DT,n>& x, const Entity& e, 
		   const FieldVector<DT,n>& xi) const = 0; 
	  
	virtual RT S0 (const FieldVector<DT,n>& x, const Entity& e, 
			const FieldVector<DT,n>& xi) const
	{}
	
	virtual RT C1_0 (const FieldVector<DT,n>& x, const Entity& e, 
			const FieldVector<DT,n>& xi) const = 0;
	  
	virtual const FieldVector<DT,n>& vTotal (const Entity& e, const int numberInSelf) = 0;
	
	virtual RT press (const FieldVector<DT,n>& x, const Entity& e, 
			const FieldVector<DT,n>& xi) const = 0;
	
	virtual RT pressBC (const FieldVector<DT,n>& x, const Entity& e, 
			const FieldVector<DT,n>& xi) const = 0;
	
	virtual RT porosity ()
	{
		return 1.0;
	}

	TransportProblem2p2c(Henry& h, TwoPhaseRelations& law = *(new LinearLaw), 
								 const bool cap = false) 
	: capillary(cap), materialLaw(law), henry(h)
	{	
	}
	
	virtual ~TransportProblem2p2c () {}
	
	Henry& henry;
	const bool capillary;
	TwoPhaseRelations& materialLaw;
	VelType velocity;
  };

}
#endif
