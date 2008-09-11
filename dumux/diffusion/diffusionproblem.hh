// $Id$ 

#ifndef DUNE_DIFFUSIONPROBLEM_HH
#define DUNE_DIFFUSIONPROBLEM_HH

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
 * @brief  Base class for defining an instance of the diffusion problem
 * @author Bernd Flemisch
 */
/**
 * \ingroup diffusion
 * \defgroup diffusionProblems Problems
 */

namespace Dune
{
  /*! \ingroup diffusionProblems 
   * @brief base class that defines the parameters of a diffusion equation 
   * 
   * An interface for defining parameters for the stationary diffusion equation 
   * \f$ - \text{div}\, (\lambda K \text{grad}\, p ) = q, \f$, 
   * \f$p = g\f$ on \f$\Gamma_1\f$, and \f$\lambda K \text{grad}\, p = J\f$ 
   * on \f$\Gamma_2\f$. Here, 
   * \f$p\f$ denotes the pressure, \f$K\f$ the absolute permeability, 
   * and \f$\lambda\f$ the total mobility, possibly depending on the 
   * saturation. 
   *
   *	Template parameters are:
   *	
   *	- Grid  a DUNE grid type
   *	- RT    type used for return values 
   */
  template<class G, class RT, class VC>
  class DiffusionProblem {
  protected:
	typedef typename G::ctype DT;
	typedef typename G::Traits::template Codim<0>::Entity Entity;

  public:
	  typedef G GridType;
	  typedef RT ReturnType;
	  enum {n=G::dimension, m=1};
	  
	//! evaluate diffusion tensor
	/*! Evaluate the diffusion tensor at given location
	  @param[in]  x    position in global coordinates
	  @param[in]  e    entity of codim 0
	  @param[in]  xi   position in reference element of e
	  @param[out] D    diffusion tensor to be filled
	 */
	virtual FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) = 0;
	
	//! evaluate saturation
	/*! Evaluate the saturation at given location
	  @param[in]  x    position in global coordinates
	  @param[in]  e    entity of codim 0
	  @param[in]  xi   position in reference element of e
	 */

	//! evaluate source term
	/*! evaluate source term at given location
	  @param[in]  x    position in global coordinates
	  @param[in]  e    entity of codim 0
	  @param[in]  xi   position in reference element of e
	  \return     value of source term
	 */
	virtual RT q   (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) = 0;

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
	
	virtual RT gSat (const FieldVector<DT,n>& x, const Entity& e, 
					  const FieldVector<DT,n>& xi) const 
	{
		return 1;
	}
	  
	//! evaluate Neumann boundary condition at given position
	/*! evaluate Neumann boundary condition at given position
	  @param[in]  x    position in global coordinates
	  \return     boundary condition value
	 */
	virtual RT J (const FieldVector<DT,n>& x, const Entity& e, 
				  const FieldVector<DT,n>& xi) const = 0;
	  
	const FieldVector<DT,n>& gravity() const
	{
		return gravity_;
	}
	
	//! constructor
	/** @param law implementation of Material laws. Class TwoPhaseRelations or derived.
	 *  @param cap flag to include capillary forces. 
	 */
	DiffusionProblem(VC& variableobject, TwoPhaseRelations& law = *(new LinearLaw), 
			const bool cap = false, FieldVector<DT,n> g = *(new FieldVector<DT,n>(0))) 
	: variables(variableobject), materialLaw(law), capillary(cap), gravity_(g)
	{	}
	
	//! always define virtual destructor in abstract base class
	virtual ~DiffusionProblem () {}
	
	//! a class describing relations between two phases and the porous medium
	VC& variables;
	TwoPhaseRelations& materialLaw;
	const bool capillary;
  protected:
	  FieldVector<DT,n> gravity_;
  };

}
#endif
