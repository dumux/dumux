// $Id$ 

#ifndef DUNE_STOKESPROBLEM_HH
#define DUNE_STOKESPROBLEM_HH

#include<iostream>
#include<iomanip>

#include<dune/common/fvector.hh>
#include<dune/common/fmatrix.hh>
#include<dune/istl/bvector.hh>
#include<dune/common/exceptions.hh>
#include<dune/grid/common/grid.hh>
#include<dune/grid/common/referenceelements.hh>
#include<dune/grid/utility/intersectiongetter.hh>
#include<dune/disc/operators/boundaryconditions.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the Stokes problem
 * @author Bernd Flemisch
 */

namespace Dune {
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
template<class G, class RT> class StokesProblem {
	typedef typename G::ctype DT;
	enum {dim=G::dimension, m=G::dimension+1};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator
			IntersectionIterator;

public:
	//! evaluate source term of the momentum equation 
	/*! evaluate source term of the momentum equation at given location
	 @param[in]  x    position in global coordinates
	 @param[in]  e    entity of codim 0
	 @param[in]  xi   position in reference element of e
	 \return     value of source term
	 */
	virtual FieldVector<RT,dim> q(const FieldVector<DT,dim>& x, const Entity& e,
			const FieldVector<DT,dim>& xi) const = 0;

	//! return type of boundary condition at the given global coordinate
	/*! return type of boundary condition at the given global coordinate
	 @param[in]  x    position in global coordinates
	 \return     boundary condition type given by enum in this class
	 */
	virtual BoundaryConditions::Flags bctype (const FieldVector<DT,dim>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,dim>& xi) const = 0;

	//! evaluate velocity Dirichlet boundary condition at given position
	/*! evaluate velocity Dirichlet boundary condition at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual FieldVector<RT,dim> g(const FieldVector<DT,dim>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,dim>& xi) const = 0;

	//! evaluate normal force boundary condition at given position
	/*! evaluate normal force boundary condition \f$ p - \mu \vec{n}\cdot(\nabla\vec{u}\cdot\vec{n}) = J_n \f$ at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual RT Jn(const FieldVector<DT,dim>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,dim>& xi) const
        {
	  DUNE_THROW(NotImplemented, "no Jn specified, but requested");
	  
	  return 0; 
	}

	//! evaluate tangential force boundary condition at given position
	/*! evaluate tangential force boundary condition \f$- \mu (\nabla\vec{u}\cdot\vec{n})_\tau = \vec{J}_\tau \f$ at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual FieldVector<RT,dim> Jt(const FieldVector<DT,dim>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,dim>& xi) const
        {
	  DUNE_THROW(NotImplemented, "no Jt specified, but requested");
	  
	  FieldVector<RT,dim> result(0);
	  return result; 
	}

	//! evaluate Beavers-Joseph proportionality constant at given position
	/*! evaluate Beavers-Joseph proportionality constant \f$c = \sqrt(k)/\alpha\f$ 
	  such that \f$u_\tau = - c (\nabla u\cdot n)_\tau\f$ 
	 @param[in]  x    position in global coordinates
	 \return     value of the proportionality constant 
	 */
	virtual RT beaversJosephC(const FieldVector<DT,dim>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,dim>& xi) const
        {
	  return 0;
	}

	//! evaluate viscosity at given position
	/*! evaluate viscosity at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual RT mu(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const = 0;

  virtual FieldVector<RT,dim> velocity(const FieldVector<DT,dim>& x) const
  {
    DUNE_THROW(NotImplemented, "no exact solution available");

    FieldVector<RT,dim> result(0);
    return result; 
  }

  virtual RT pressure(const FieldVector<DT,dim>& x) const
  {
    DUNE_THROW(NotImplemented, "no exact solution available");

    return 0;
  }

  virtual FieldMatrix<DT, dim, dim> velocityGradient(const FieldVector<DT,dim>& x) const
  {
    DUNE_THROW(NotImplemented, "no exact solution available");

    FieldMatrix<DT, dim, dim> result(0); 
    return result;
  }

	StokesProblem()
        {}

	//! always define virtual destructor in abstract base class
	virtual ~StokesProblem() 
        {}

protected:
};

}
#endif
