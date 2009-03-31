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
 *    Template parameters are:
 *
 *    - Grid  a DUNE grid type
 *    - Scalar    type used for return values
 */
template<class Grid, class Scalar> class StokesProblem
{
    enum {dim=Grid::dimension, numEq=Grid::dimension+1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::template Codim<0>::LeafIntersectionIterator::IntersectionIterator
    IntersectionIterator;

public:
    //! evaluate source term of the momentum equation
    /*! evaluate source term of the momentum equation at given location
      @param[in]  globalPos    position in global coordinates
      @param[in]  element    entity of codim 0
      @param[in]  localPos   position in reference element
      \return     value of source term
    */
    virtual FieldVector<Scalar,numEq> q(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                        const FieldVector<Scalar,dim>& localPos) const
    {
        DUNE_THROW(NotImplemented, "no q specified, but requested");

        FieldVector<Scalar,numEq> result(0);
        return result;
    }

    //! return type of boundary condition at the given global coordinate
    /*! return type of boundary condition at the given global coordinate
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition type given by enum in this class
    */
    virtual BoundaryConditions::Flags bctype (const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                              const IntersectionIterator& intersectionIt,
                                              const FieldVector<Scalar,dim>& localPos) const = 0;

    //! evaluate velocity Dirichlet boundary condition at given position
    /*! evaluate velocity Dirichlet boundary condition at given position
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition value
    */
    virtual FieldVector<Scalar,numEq> g(const FieldVector<Scalar,dim>& globalPos, const Element& element,
					const IntersectionIterator& intersectionIt,
					const FieldVector<Scalar,dim>& localPos) const = 0;

    //! evaluate Neumann boundary condition at given position
    /*! evaluate Neumann boundary condition at given position
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition value
    */
    virtual FieldVector<Scalar,numEq> J(const FieldVector<Scalar,dim>& globalPos, const Element& element,
					const IntersectionIterator& intersectionIt,
					const FieldVector<Scalar,dim>& localPos)
    {
        DUNE_THROW(NotImplemented, "no J specified, but requested");

        FieldVector<Scalar,numEq> result(0);
        return result;
    }

    //! evaluate normal force boundary condition at given position
    /*! evaluate normal force boundary condition \f$ p - \mu \vec{n}\cdot(\nabla\vec{u}\cdot\vec{n}) = J_n \f$ at given position
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition value
    */
    virtual Scalar Jn(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                      const IntersectionIterator& intersectionIt,
                      const FieldVector<Scalar,dim>& localPos) const
    {
        DUNE_THROW(NotImplemented, "no Jn specified, but requested");

        return 0;
    }

    //! evaluate tangential force boundary condition at given position
    /*! evaluate tangential force boundary condition \f$- \mu (\nabla\vec{u}\cdot\vec{n})_\tau = \vec{J}_\tau \f$ at given position
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition value
    */
    virtual FieldVector<Scalar,dim> Jt(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                       const IntersectionIterator& intersectionIt,
                                       const FieldVector<Scalar,dim>& localPos) const
    {
        DUNE_THROW(NotImplemented, "no Jt specified, but requested");

        FieldVector<Scalar,dim> result(0);
        return result;
    }

    //! evaluate Beavers-Joseph proportionality constant at given position
    /*! evaluate Beavers-Joseph proportionality constant \f$c = \sqrt(k)/\alpha\f$
      such that \f$u_\tau = - c (\mu \nabla u\cdot n)_\tau\f$
      @param[in]  globalPos    position in global coordinates
      \return     value of the proportionality constant
    */
    virtual Scalar beaversJosephC(const FieldVector<Scalar,dim>& globalPos, const Element& element,
                                  const IntersectionIterator& intersectionIt,
                                  const FieldVector<Scalar,dim>& localPos) const
    {
        return 0;
    }

    //! evaluate viscosity at given position
    /*! evaluate viscosity at given position
      @param[in]  globalPos    position in global coordinates
      \return     boundary condition value
    */
    virtual Scalar mu(const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos) const = 0;

    virtual FieldVector<Scalar,numEq> velocity(const FieldVector<Scalar,dim>& globalPos) const
    {
        DUNE_THROW(NotImplemented, "no exact solution available");

        FieldVector<Scalar,numEq> result(0);
        return result;
    }

    virtual Scalar pressure(const FieldVector<Scalar,dim>& globalPos) const
    {
        DUNE_THROW(NotImplemented, "no exact solution available");

        return 0;
    }

    virtual FieldMatrix<Scalar, dim, dim> velocityGradient(const FieldVector<Scalar,dim>& globalPos) const
    {
        DUNE_THROW(NotImplemented, "no exact solution available");

        FieldMatrix<Scalar, dim, dim> result(0);
        return result;
    }

    virtual void dirichletIndex(const Dune::FieldVector<Scalar,dim>& globalPos, const Element& element,
                                const IntersectionIterator& intersectionIt,
                                const Dune::FieldVector<Scalar,dim>& localPos, Dune::FieldVector<int,numEq>& dirichletIndex) const
    {
        for (int i = 0; i < numEq; i++)
            dirichletIndex[i]=i;
        return;
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
