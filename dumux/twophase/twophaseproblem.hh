// $Id$ 

#ifndef DUNE_TWOPHASEPROBLEM_HH
#define DUNE_TWOPHASEPROBLEM_HH

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
#include<dumux/material/twophaserelations_deprecated.hh>
#include<dumux/material/linearlaw_deprecated.hh>

/**
 * @file
 * @brief  Base class for defining an instance of the TwoPhase problem
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
template<class G, class RT> class TwoPhaseProblem {
	typedef typename G::ctype DT;
	enum {n=G::dimension, m=2};
	typedef typename G::Traits::template Codim<0>::Entity Entity;
	typedef typename IntersectionIteratorGetter<G,LeafTag>::IntersectionIterator
			IntersectionIterator;

public:
	//! evaluate diffusion tensor
	/*! Evaluate the diffusion tensor at given location
	 @param[in]  x    position in global coordinates
	 @param[in]  e    entity of codim 0
	 @param[in]  xi   position in reference element of e
	 @param[out] D    diffusion tensor to be filled
	 */
	virtual const FieldMatrix<DT,n,n>& K(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) = 0;

	//! evaluate source term
	/*! evaluate source term at given location
	 @param[in]  x    position in global coordinates
	 @param[in]  e    entity of codim 0
	 @param[in]  xi   position in reference element of e
	 \return     value of source term
	 */
	virtual FieldVector<RT,m> q(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const = 0;

	//! return type of boundary condition at the given global coordinate
	/*! return type of boundary condition at the given global coordinate
	 @param[in]  x    position in global coordinates
	 \return     boundary condition type given by enum in this class
	 */
	//	virtual FieldVector<BoundaryConditions::Flags, m> bctype (const FieldVector<DT,n>& x, const Entity& e, 
	//			const IntersectionIterator& intersectionIt, 
	//			const FieldVector<DT,n>& xi) const = 0;

	virtual FieldVector<BoundaryConditions::Flags, m>bctype(
			const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi) const = 0;

	//! returns index of the primary variable corresponding to the dirichlet boundary condition at the given global coordinate
		/*! returns index of the primary variable corresponding to the dirichlet boundary condition at the given global coordinate
		 @param[in]  x    position in global coordinates
		 \return     index of the primary variable
		 */
	
	virtual void dirichletIndex(const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi, FieldVector<int,m>& dirichletIdx) const 
	{
		for (int i = 0; i < m; i++)
			dirichletIdx[i]=i;
		return;
	}

	//! evaluate Dirichlet boundary condition at given position
	/*! evaluate Dirichlet boundary condition at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual FieldVector<RT,m> g(const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi) const = 0;

	//! evaluate Neumann boundary condition at given position
	/*! evaluate Neumann boundary condition at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual FieldVector<RT,m> J(const FieldVector<DT,n>& x, const Entity& e,
			const IntersectionIterator& intersectionIt,
			const FieldVector<DT,n>& xi) const = 0;

	//! evaluate initial condition at given position
	/*! evaluate initial boundary condition at given position
	 @param[in]  x    position in global coordinates
	 \return     boundary condition value
	 */
	virtual FieldVector<RT,m> initial(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) const = 0;

	virtual double porosity(const FieldVector<DT,n>& x, const Entity& e,
			const FieldVector<DT,n>& xi) const = 0;

	virtual FieldVector<RT,n> gravity() const = 0;

	virtual FieldVector<RT,4> materialLawParameters(const FieldVector<DT,n>& x,
			const Entity& e, const FieldVector<DT,n>& xi) const = 0;

	TwoPhaseRelations& materialLaw() {
		return materialLaw_;
	}

	//element-wise return of the values of an Exact solution
	virtual RT uExOutVertex(int &ElementIndex, int VariableIndex) const {
		DUNE_THROW(NotImplemented, "Ex(akt) Solution");
		return 0;
	}
	
	//updates an exact/analytic solution
	virtual void updateExSol(double &dt,
			BlockVector<FieldVector<RT, m> > &approxSol) {
		DUNE_THROW(NotImplemented, "Ex(akt) Solution");
		return;
	}

	TwoPhaseProblem(TwoPhaseRelations& law = *(new
			LinearLaw), const bool exsol = false) :
	  exsolution(exsol),  materialLaw_(law) {
	}

	//! always define virtual destructor in abstract base class
	virtual ~TwoPhaseProblem() {
	}

	const bool exsolution;

protected:
	TwoPhaseRelations& materialLaw_;
};

}
#endif
