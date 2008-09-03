#ifndef DUNE_FRACTIONALFLOWPROBLEM_HH
#define DUNE_FRACTIONALFLOWPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations.hh>

#include <dumux/material/property_baseclasses.hh>
#include "dumux/fractionalflow/variableclass.hh"

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
  class FractionalFlowProblem {
  protected:
	typedef typename G::ctype DT;
	typedef typename G::Traits::template Codim<0>::Entity Entity;

  public:
	  typedef G GridType;
	  typedef RT ReturnType;
	  enum {n=G::dimension, m=1};

	//! evaluate source term for the pressure equation
	/*! evaluate source term for the pressure equation at given location
	  @param[in]  x    position in global coordinates
	  @param[in]  e    entity of codim 0
	  @param[in]  xi   position in reference element of e
	  \return     value of source term
	 */
	virtual RT qPress  (const FieldVector<DT,n>& x, const Entity& e,
					const FieldVector<DT,n>& xi) = 0;

	//! return type of boundary condition for the pressure equation at the given global coordinate
	/*! return type of boundary condition for the pressure equation at the given global coordinate
	  @param[in]  x    position in global coordinates
	  \return     boundary condition type given by enum in this class
	 */
	virtual BoundaryConditions::Flags bctypePress (const FieldVector<DT,n>& x, const Entity& e,
					   const FieldVector<DT,n>& xi) const = 0;

	//! return type of boundary condition for the saturation equation at the given global coordinate
	/*! return type of boundary condition for the saturation equation at the given global coordinate
	  @param[in]  x    position in global coordinates
	  \return     boundary condition type given by enum in this class
	 */
	virtual BoundaryConditions::Flags bctypeSat (const FieldVector<DT,n>& x, const Entity& e,
					   const FieldVector<DT,n>& xi) const = 0;

	//! evaluate Dirichlet boundary condition for the pressure equation at given position
	/*! evaluate Dirichlet boundary condition for the pressure equation at given position
	  @param[in]  x    position in global coordinates
	  \return     boundary condition value
	 */
	virtual RT gPress (const FieldVector<DT,n>& x, const Entity& e,
				  const FieldVector<DT,n>& xi) const = 0;

	//! evaluate Dirichlet boundary condition for the saturation equation at given position
	/*! evaluate Dirichlet boundary condition for the saturation equation at given position
	  @param[in]  x    position in global coordinates
	  \return     boundary condition value
	 */
	virtual RT gSat (const FieldVector<DT,n>& x, const Entity& e,
					  const FieldVector<DT,n>& xi) const
	{
		return 1;
	}

	//! evaluate Neumann boundary condition for the pressure equation at given position
	/*! evaluate Neumann boundary condition for the pressure equation at given position
	  @param[in]  x    position in global coordinates
	  \return     boundary condition value
	 */
	virtual RT JPress (const FieldVector<DT,n>& x, const Entity& e,
				  const FieldVector<DT,n>& xi) const = 0;

	//! evaluate initial condition for saturation at given position
	/*! evaluate initial condition for saturation at given position
	  @param[in]  x    position in global coordinates
	  \return    initial condition value
	 */
	virtual RT S0 (const FieldVector<DT,n>& x, const Entity& e,
				  const FieldVector<DT,n>& xi) const = 0;

	//! evaluate gravity
	/*! evaluate gravity
	  \return     gravity vector
	 */
	const FieldVector<DT,n>& gravity() const
	{
		return gravity_;
	}

	//! constructor
	/** @param law implementation of Material laws. Class TwoPhaseRelations or derived.
	 *  @param cap flag to include capillary forces.
	 */
	FractionalFlowProblem(VC& variableobject, MediumNonIsothermal& wp, MediumNonIsothermal& nwp, Matrix2p<G, RT>& s, TwoPhaseRelations<G, RT>& law = *(new TwoPhaseRelations<G,RT>), const bool cap = false)
	 : variables(variableobject), wettingphase(wp), nonwettingphase(nwp), soil(s), capillary(cap), materialLaw(law)
	 {}

	//! always define virtual destructor in abstract base class
	virtual ~FractionalFlowProblem () {}

	//! a class describing relations between two phases and the porous medium
	VC& variables;
	MediumNonIsothermal& wettingphase;
	MediumNonIsothermal& nonwettingphase;
	Matrix2p<G, RT>& soil;
	TwoPhaseRelations<G, RT>& materialLaw;
	const bool capillary;


  protected:
	  FieldVector<DT,n> gravity_;
  };

}
#endif
