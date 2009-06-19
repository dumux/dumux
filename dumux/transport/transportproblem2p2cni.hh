// $Id:$
/*****************************************************************************
* Copyright (C) 2009 by Jochen Fritz                                         *
* Institute of Hydraulic Engineering                                         *
* University of Stuttgart, Germany                                           *
* email: <givenname>.<name>@iws.uni-stuttgart.de                             *
*                                                                            *
* This program is free software; you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation; either version 2 of the License, or          *
* (at your option) any later version, as long as this copyright notice       *
* is included in its original form.                                          *
*                                                                            *
* This program is distributed WITHOUT ANY WARRANTY.                          *
*****************************************************************************/

#ifndef TRANSPORTPROBLEM2P2CNI_HH
#define TRANSPORTPROBLEM2P2CNI_HH

#include <iostream>
#include <iomanip>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/disc/operators/boundaryconditions.hh>

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/twophaserelations.hh>
#include <dumux/operators/boundaryconditions2p2c.hh>
#include <dumux/fractionalflow/variableclass2p2cni.hh>

namespace Dune
{
//! Base class for the definition of nonisothermal 2p2c problems
/** This base class defines all boundary and initial functions which are needed
 * for a decoupled nonisothermal 2p2c computation.
 */
template<class G, class Scalar>
class TransportProblem2p2cni
{
    enum {dim = G::dimension};
    typedef typename G::ctype ct;
    typedef typename G::Traits::template Codim<0>::Entity Entity;

public:
    //! Type of concentration boundary condition.
    /**    either the concentration or the saturation have to be defined
     * on boundaries with dirichlet pressure BCs.
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual BoundaryConditions2p2c::Flags bc_type (const FieldVector<ct,dim>& x, const Entity& e,
                                                   const FieldVector<ct,dim>& xi) const = 0;

    //! Type of concentration initisl condition.
    /**    either the concentration or the saturation have to be defined
     * as initial condition.
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual BoundaryConditions2p2c::Flags initcond_type (const FieldVector<ct,dim>& x, const Entity& e,
                                                  const FieldVector<ct,dim>& xi) const = 0;

    //! Type of pressure boundary condition.
    /**    Pressure (dirichlet) or flux (neumann) have to be defined on boundaries.
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual BoundaryConditions::Flags press_bc_type (const FieldVector<ct,dim>& x, const Entity& e,
                                               const FieldVector<ct,dim>& xi) const = 0;

    virtual BoundaryConditions::Flags heat_bc_type (const FieldVector<ct,dim>& x, const Entity& e,
                                               const FieldVector<ct,dim>& xi) const = 0;

    //! Feed concentration boundary condition
    /** Feed concentration is the (global) mass fraction of component 1 in the mixture
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual Scalar dirichletConcentration (const FieldVector<ct,dim>& x, const Entity& e,
                   const FieldVector<ct,dim>& xi) const = 0;

    //! Saturation boundary condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual Scalar dirichletSat (const FieldVector<ct,dim>& x, const Entity& e,
                   const FieldVector<ct,dim>& xi) const = 0;

    //! Pressure (dirichlet) boundary condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual Scalar dirichlet (const FieldVector<ct,dim>& x, const Entity& e,
                       const FieldVector<ct,dim>& xi) const = 0;

    //! Temperature boundary condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual Scalar dirichletTemp (const FieldVector<ct,dim>& x, const Entity& e,
                   const FieldVector<ct,dim>& xi) const = 0;

    //! Flux (neumann) boundary condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual FieldVector<Scalar,3> neumann(const FieldVector<ct,dim>& x, const Entity& e,
                                 const FieldVector<ct,dim>& xi) const = 0;

    //! Source of components
    /** Describes the source of the components per unit volume in \f$ \left[ frac{kg}{m^3} \right] \f$
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual FieldVector<Scalar,2> source (const FieldVector<ct,dim>& x, const Entity& e,
                                 const FieldVector<ct,dim>& xi) const = 0;

    //! Heat source
    /** Heat flux from outside per unit volume in \f$ \left[ frac{W}{m^3} \right] \f$
     * @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual Scalar source_heat (const FieldVector<ct,dim>& x, const Entity& e,
                   const FieldVector<ct,dim>& xi) const = 0;

    //! Saturation initial condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual Scalar initSat (const FieldVector<ct,dim>& x, const Entity& e,
                   const FieldVector<ct,dim>& xi) const = 0;

    //! Feed concentration initial condition
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual Scalar initConcentration (const FieldVector<ct,dim>& x, const Entity& e,
                     const FieldVector<ct,dim>& xi) const = 0;

    //! Temperature initial condition \f$ \left[ K \right] \f$
    /** @param x global coordinates
     * @param e reference to the cell for which the function is to be evaluated
     * @param xi local coordinates inside e
     */
    virtual Scalar initTemp (const FieldVector<ct,dim>& x, const Entity& e,
                    const FieldVector<ct,dim>& xi) const = 0;

    //! gravity vector
    virtual const FieldVector<Scalar,dim> gravity()
    {
        FieldVector<Scalar,dim> gravity_(0.);
        return gravity_;
    }


    //! Constructor
    /**
     *
     */
    TransportProblem2p2cni(Dune::VariableClass2p2cni<G, Scalar>& var, Liquid_GL& liq, Gas_GL& gas, Matrix2p<G, Scalar>& s, TwoPhaseRelations<G, Scalar>& law)
        :variables(var), liquidPhase(liq), gasPhase(gas), soil(s), materialLaw(law)
    {
    }

    virtual ~TransportProblem2p2cni ()
    {
    }

    TwoPhaseRelations<G, Scalar>& materialLaw;
    Liquid_GL& liquidPhase;
    Gas_GL& gasPhase;
    Matrix2p<G, Scalar>& soil;
    VariableClass2p2cni<G, Scalar>& variables;
};

}
#endif
