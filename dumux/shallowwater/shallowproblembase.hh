// $Id: transportproblem.hh 689 2008-10-10 15:48:03Z lauser $

#ifndef DUNE_SHALLOWPROBLEMBASE_HH
#define DUNE_SHALLOWPROBLEMBASE_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/shallowwater/solidsurfacebase.hh>
#include<dumux/shallowwater/shallowvariableclass.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the transport problem
 * @author Bernd Flemisch
 */

namespace Dune
{
/**  \ingroup transport
 *   \defgroup transportProblems Problems
 */
/*! \ingroup transportProblems
 *   *- Template parameters are:
 *   *- Grid  a DUNE grid type
 *	- DT    type used for return values
 */

template<class G, class DT, class VC> class ShallowProblemBase
{
	enum
	{	dim=G::dimension, m=1, blocksize=2*G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef Dune::FieldVector<DT, dim> VelType;
	typedef Dune::FieldVector<DT,dim+1> SystemType;


public:
	typedef Dune::SolidSurfaceBase<G,DT> Surface;
	
	/*! return type of boundary condition at the given global coordinate
	 @param[in]  x    position in global coordinates
	 \return     boundary condition type given by enum in this class
	 */
	virtual BoundaryConditions::Flags bctype(const FieldVector<DT,dim>& x,
			const Entity& e, const FieldVector<DT,dim>& xi) const = 0;

	/*! evaluate Dirichlet boundary condition at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual DT dirichlet(const FieldVector<DT,dim>& x, const Entity& e,
			const FieldVector<DT,dim>& xi)const = 0;

	/*! evaluate Neumann boundary condition at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual SystemType neumann(const FieldVector<DT,dim>& x, const Entity& e,
			const FieldVector<DT,dim>& xi) const = 0;
	
	/*! evaluate initial boundary condition at given position
	 @param[in]  x position in global coordinates
	 \return    initial condition value
	 */
	virtual DT setInitWDepth(const FieldVector<DT,dim>& x, const Entity& e,
			const FieldVector<DT,dim>& xi) const = 0;

	virtual VelType setInitVel(const FieldVector<DT,dim>& x, const Entity& e,
			const FieldVector<DT,dim>& xi) const = 0; 

	//Definition of sources

	virtual DT setSource(const FieldVector<DT,dim>& x, const Entity& e,
			const FieldVector<DT,dim>& xi) const = 0;

	ShallowProblemBase(VC& variableobject, Surface& surfaceobject) :
		variables(variableobject), surface(surfaceobject)
	{
	}

	//! always define virtual destructor in abstract base class
	virtual ~ShallowProblemBase()
	{
	}

	VC& variables;
	Surface& surface;

};
}
#endif
