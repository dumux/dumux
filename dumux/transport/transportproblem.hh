#ifndef DUNE_TRANSPORTPROBLEM_HH
#define DUNE_TRANSPORTPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations_deprecated.hh>
#include<dumux/material/linearlaw_deprecated.hh>

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
 //! \ingroup transportProblems
 /**  @brief  base class that defines the parameters of a transport equation
   *   An interface for defining parameters for the scalar transport equation 
   *  \f$S_t - \text{div}\, (f_\text{w}(S) \boldsymbol{v}_\text{total}) = 0\f$, 
   * \f$S = g\f$ on \f$\Gamma_1\f$, and \f$S(t = 0) = S_0\f$. Here, 
   * \f$S\f$ denotes the wetting phase saturation, 
   * \f$\boldsymbol{v}_\text{total}\f$ the total velocity, 
   * and \f$f_\text{w}\f$ the wetting phase fractional flow function.
   *
   *	Template parameters are:
   *	
   *	- Grid  a DUNE grid type
   *	- RT    type used for return values 
   */
  template<class G, class RT, class VC>
  class TransportProblem {
	typedef typename G::ctype DT;
	enum {n=G::dimension, m=1, blocksize=2*G::dimension};
	typedef typename G::Traits::template Codim<0>::Entity Entity;

  public:
	//! return type of boundary condition at the given global coordinate
	/*! return type of boundary condition at the given global coordinate
	  @param[in]  x    position in global coordinates
	  \return     boundary condition type given by enum in this class
	 */
	virtual BoundaryConditions::Flags bctype (const FieldVector<DT,n>& x, const Entity& e, 
					   const FieldVector<DT,n>& xi) const = 0;

	//! evaluate Dirichlet boundary condition at given position
	/*! evaluate Dirichlet boundary condition at given position
	  @param[in]  x    position in global coordinates
	  \return     boundary condition value
	 */
	virtual RT g (const FieldVector<DT,n>& x, const Entity& e, 
				  const FieldVector<DT,n>& xi) const = 0;
	  
	//! evaluate initial condition at given position
	/*! evaluate initial boundary condition at given position
	  @param[in]  x    position in global coordinates
	  \return    initial condition value
	 */
	virtual RT S0 (const FieldVector<DT,n>& x, const Entity& e, 
				  const FieldVector<DT,n>& xi) const = 0;
	
	const FieldVector<DT,n>& gravity()
	{	  
		FieldVector<DT,n> gravity_ = 0;
		return gravity_;
	}
	
//	virtual const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
//					const FieldVector<DT,n>& xi) = 0;
	  
	//! evaluate velocity 
	/*! Evaluate the velocity at the element faces
	  @param[in]  e              entity of codim 0
	  @param[in]  numberInSelf   local index of element face
	  @param[out] vTotal         velocity vector to be filled
	 */


	//! Returns the porosity of the porous medium
	virtual RT porosity () const
	{
		return 1.0;
	}
	
	virtual BlockVector<FieldVector<RT, 2> >& getuEx()
	{
		DUNE_THROW(NotImplemented, "Ex(akt) Solution");
		return uE;
	}
	
	//updates an exact/analytic solution
	virtual void updateExSol() {
		DUNE_THROW(NotImplemented, "Ex(akt) Solution");
		return;
	}

	virtual void settime(double &dt)
	{
		DUNE_THROW(NotImplemented, "Ex(akt) Solution");
		return;	
	}
	
	//! constructor
	/** @param law implementation of material laws. Class TwoPhaseRelations or derived.
	 *  @param cap flag for including capillary forces.
	 */

	TransportProblem(VC& variableobject, TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false, const bool exsol = false) 
	: variables(variableobject), capillary(cap), materialLaw(law),exsolution(exsol), uE(0)
	{	}
	
	//! always define virtual destructor in abstract base class
	virtual ~TransportProblem () {}
	
	const bool exsolution;
	const bool capillary;
	TwoPhaseRelations& materialLaw;
	VC& variables;
	BlockVector<FieldVector<RT, 2> > uE;
  };

}
#endif
