#ifndef DUNE_PARABOLICPROBLEM_HH
#define DUNE_PARABOLICPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/disc/operators/boundaryconditions.hh>
#include<dumux/material/twophaserelations.hh>
#include<dumux/material/linearlaw.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the Parabolic problem
 * @author Bernd Flemisch
 */

namespace Dune
{
  //! base class that defines the parameters of a diffusion equation
  /*! An interface for defining parameters for the stationary diffusion equation 
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
  template<class G, class RT>
  class ParabolicProblem {
	typedef typename G::ctype DT;
	enum {n=G::dimension, m=1};
	typedef typename G::Traits::template Codim<0>::Entity Entity;

  public:
	//! evaluate diffusion tensor
	/*! Evaluate the diffusion tensor at given location
	  @param[in]  x    position in global coordinates
	  @param[in]  e    entity of codim 0
	  @param[in]  xi   position in reference element of e
	  @param[out] D    diffusion tensor to be filled
	 */
	virtual const FieldMatrix<DT,n,n>& K (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) = 0;

	//! evaluate source term
	/*! evaluate source term at given location
	  @param[in]  x    position in global coordinates
	  @param[in]  e    entity of codim 0
	  @param[in]  xi   position in reference element of e
	  \return     value of source term
	 */
	virtual RT q   (const FieldVector<DT,n>& x, const Entity& e, 
					const FieldVector<DT,n>& xi) const = 0;

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
	  
	//! evaluate Neumann boundary condition at given position
	/*! evaluate Neumann boundary condition at given position
	  @param[in]  x    position in global coordinates
	  \return     boundary condition value
	 */
	virtual RT J (const FieldVector<DT,n>& x, const Entity& e, 
				  const FieldVector<DT,n>& xi) const = 0;
	  
	//! evaluate initial condition at given position
	/*! evaluate initial boundary condition at given position
	  @param[in]  x    position in global coordinates
	  \return     boundary condition value
	 */
	virtual RT initial (const FieldVector<DT,n>& x, const Entity& e, 
				  const FieldVector<DT,n>& xi) const = 0;
	  
	ParabolicProblem(TwoPhaseRelations& law = *(new LinearLaw), const bool cap = false) 
	: materialLaw(law), capillary(cap)
	{	}
	
	//! always define virtual destructor in abstract base class
	virtual ~ParabolicProblem () {}
	
	TwoPhaseRelations& materialLaw;
	const bool capillary;
  };

}
#endif
