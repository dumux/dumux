// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch, Jochen Fritz, Markus Wolff   *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/

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

/**
 * @file
 * @brief  Base class for defining an instance of the diffusion problem
 * @author Bernd Flemisch
 */

namespace Dune
{
/*! \ingroup fracflow
 * @brief base class that defines the parameters of loosely coupled diffusion and transport equations
 *
 * An interface for defining parameters for the stationary diffusion equation
 *  \f[\text{div}\, \boldsymbol{v} = q\f]
 *  and a scalar transport equation
 *  \f[
 *    \frac{\partial S}{\partial t} + \text{div}\, \boldsymbol{v_\alpha} = 0,
 *  \f]
 *  where, the velocity \f$\boldsymbol{v} \sim \boldsymbol{K} \nabla p \f$,
 *  \f$p\f$ is a pressure and q a source/sink term, \f$S\f$ denotes a phase saturation and \f$\boldsymbol{v_\alpha}\f$ is a phase velocity.
 *
 *  Template parameters are:
 *
 *  - GridView      a DUNE gridview type
 *  - Scalar        type used for scalar quantities
 *  - VC            type of a class containing different variables of the model
 */
template<class GridView, class Scalar, class VC>
class FractionalFlowProblem
{
public:
    enum
        {    dim=GridView::dimension,dimWorld= GridView::dimensionworld,numEq=1};
protected:
    typedef typename GridView::Grid Grid;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    //! evaluate source term
    /*! evaluate source term at given location
     @param  globalPos    position in global coordinates
     @param  element      entity of codim 0
     @param  localPos     position in reference element of element
     \return     value of source term
     */
    virtual Scalar sourcePress (const GlobalPosition& globalPos, const Element& element,
                                const LocalPosition& localPos) = 0;

    //! return type of boundary condition at the given global coordinate
    /*! return type of boundary condition at the given global coordinate
     @param  globalPos    position in global coordinates
     @param  element      entity of codim 0
     @param  localPos     position in reference element of element
     \return     boundary condition type given by enum in this class
     */
    virtual BoundaryConditions::Flags bctypePress (const GlobalPosition& globalPos, const Element& element,
                                                   const LocalPosition& localPos) const = 0;
    //! return type of boundary condition at the given global coordinate
    /*! return type of boundary condition at the given global coordinate
     @param  globalPos    position in global coordinates
     @param  element      entity of codim 0
     @param  localPos     position in reference element of element
     \return     boundary condition type given by enum in this class
     */
    virtual BoundaryConditions::Flags bctypeSat (const GlobalPosition& globalPos, const Element& element,
                                                 const LocalPosition& localPos) const = 0;

    //! evaluate Dirichlet boundary condition at given position
    /*! evaluate Dirichlet boundary condition at given position
     @param  globalPos    position in global coordinates
     @param  element      entity of codim 0
     @param  localPos     position in reference element of element
     \return     boundary condition value of a dirichlet pressure boundary condition
     */
    virtual Scalar dirichletPress (const GlobalPosition& globalPos, const Element& element,
                                   const LocalPosition& localPos) const = 0;

    //! evaluate Dirichlet boundary condition at given position
    /*! evaluate Dirichlet boundary condition at given position
     @param  globalPos    position in global coordinates
     @param  element      entity of codim 0
     @param  localPos     position in reference element of element
     \return     boundary condition value of a dirichlet saturation boundary condition
     */
    virtual Scalar dirichletSat (const GlobalPosition& globalPos, const Element& element,
                                 const LocalPosition& localPos) const
    {
        return 1;
    }

    //! evaluate Neumann boundary condition at given position
    /*! evaluate Neumann boundary condition at given position
     @param  globalPos    position in global coordinates
     @param  element      entity of codim 0
     @param  localPos     position in reference element of element
     \return     boundary condition value of a neumann pressure boundary condition.
     */
    virtual Scalar neumannPress (const GlobalPosition& globalPos, const Element& element,
                                 const LocalPosition& localPos) const = 0;

    //! evaluate Neumann boundary condition at given position
    /*! evaluate Neumann boundary condition at given position
     @param  globalPos    position in global coordinates
     @param  element      entity of codim 0
     @param  localPos     position in reference element of element
     \return     boundary condition value of a neumann saturation boundary condition.
     */
    virtual Scalar neumannSat (const GlobalPosition& globalPos, const Element& element,
                               const LocalPosition& localPos, Scalar factor) const
    {
        return 0;
    }

    //! evaluate initial condition at given position
    /*! evaluate initial boundary condition at given position
     @param  globalPos    position in global coordinates
     @param  element      entity of codim 0
     @param  localPos     position in reference element of element
     \return    initial saturation distribution.
     */
    virtual Scalar initSat (const GlobalPosition& globalPos, const Element& element,
                            const LocalPosition& localPos) const = 0;

    //! gravity constant
    /*! gravity constant
     \return    gravity vector
     */
    virtual const FieldVector<Scalar,dimWorld>& gravity() const
    {
        return gravity_;
    }

    //! properties of the soil
    /*! properties of the soil
      \return    soil
    */
    virtual Matrix2p<Grid, Scalar>& soil () const
    {
        return soil_;
    }

    //! object for definition of material law
    /*! object for definition of material law (e.g. Brooks-Corey, Van Genuchten, ...)
      \return    material law
    */
    virtual TwoPhaseRelations<Grid, Scalar>& materialLaw () const
    {
        return materialLaw_;
    }

    //! object containing different variables
    /*! object containing different variables like the primary pressure and saturation, material laws,...
     \return    variables object
     */
    virtual VC& variables ()
    {
        return variables_;
    }

    //! object containing the wetting phase parameters
    /*! object containing the wetting phase parameters (density, viscosity, ...)
     \return    wetting phase
     */
    virtual Fluid& wettingPhase () const
    {
        return wettingPhase_;
    }

    //! object containing the non-wetting phase parameters
    /*! object containing the non-wetting phase parameters (density, viscosity, ...)
     \return    non-wetting phase
     */
    virtual Fluid& nonWettingPhase () const
    {
        return nonWettingPhase_;
    }

    //! Constructs an object of type FractionalFlowProblem
    /** @param variables object of class VariableClass.
     *  @param wettingPhase implementation of a wetting phase.
     *  @param nonWettingPhase implementation of a non-wetting phase.
     *  @param soil implementation of the solid matrix
     *  @param materialLaw implementation of Material laws. Class TwoPhaseRelations or derived.
     */
    FractionalFlowProblem(VC& variables, Fluid& wettingPhase, Fluid& nonWettingPhase, Matrix2p<Grid, Scalar>& soil, TwoPhaseRelations<Grid, Scalar>& materialLaw = *(new TwoPhaseRelations<Grid,Scalar>))
        : variables_(variables), wettingPhase_(wettingPhase), nonWettingPhase_(nonWettingPhase), soil_(soil), materialLaw_(materialLaw),gravity_(0)
    {}

    //! always define virtual destructor in abstract base class
    virtual ~FractionalFlowProblem ()
    {}

private:
    VC& variables_;//object of type Dune::VariableClass
    Fluid& wettingPhase_;//object derived from Dune::Fluid
    Fluid& nonWettingPhase_;//object derived from Dune::Fluid
    Matrix2p<Grid, Scalar>& soil_;//object derived from Dune::Matrix2p
    TwoPhaseRelations<Grid, Scalar>& materialLaw_;//object of type Dune::TwoPhaseRelations or derived
    FieldVector<Scalar,dimWorld> gravity_;//vector including the gravity constant

};

}
#endif
