/*****************************************************************************
 *   Copyright (C) 2011 by Holger Class                                      *
 *   Copyright (C) 2008 by Klaus Mosthaf, Andreas Lauser, Bernd Flemisch     *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup HolleThreePThreeCModel
 */
/*!
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        holle3p3c box model.
 */
#ifndef DUMUX_HOLLE2P2C_PROPERTY_DEFAULTS_HH
#define DUMUX_HOLLE2P2C_PROPERTY_DEFAULTS_HH

#include "holle3p3cindices.hh"

#include "holle3p3cmodel.hh"
#include "holle3p3cproblem.hh"
#include "holle3p3cindices.hh"
#include "holle3p3cfluxvariables.hh"
#include "holle3p3cvolumevariables.hh"
#include "holle3p3cfluidstate.hh"
#include "holle3p3cproperties.hh"
#include "holle3p3cnewtoncontroller.hh"
// #include "holle3p3cboundaryvariables.hh"

namespace Dumux
{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(BoxHolleThreePThreeC, NumComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;

    static_assert(value == 3,
                  "Only fluid systems with 3 components are supported by the 3p3c model!");
};

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(BoxHolleThreePThreeC, NumPhases)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 3,
                  "Only fluid systems with 3 phases are supported by the 3p3c model!");
};

SET_INT_PROP(BoxHolleThreePThreeC, NumEq, 3); //!< set the number of equations to 2

//! Set the default formulation to pg-Sw-Sn
SET_INT_PROP(BoxHolleThreePThreeC,
             Formulation,
             HolleThreePThreeCFormulation::pgSwSn);

/*!
 * \brief Set the property for the material law by retrieving it from
 *        the spatial parameters.
 */
SET_PROP(BoxHolleThreePThreeC, MaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(SpatialParameters)) SpatialParameters;

public:
    typedef typename SpatialParameters::MaterialLaw type;
};

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_PROP(BoxHolleThreePThreeC, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

//! Use the holle3p3c local jacobian operator for the holle3p3c model
SET_TYPE_PROP(BoxHolleThreePThreeC,
              LocalResidual,
              HolleThreePThreeCLocalResidual<TypeTag>);

//! Use the holle3p3c specific newton controller for the holle3p3c model
SET_TYPE_PROP(BoxHolleThreePThreeC, NewtonController, HolleThreePThreeCNewtonController<TypeTag>);

//! the Model property
SET_TYPE_PROP(BoxHolleThreePThreeC, Model, HolleThreePThreeCModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(BoxHolleThreePThreeC, VolumeVariables, HolleThreePThreeCVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(BoxHolleThreePThreeC, FluxVariables, HolleThreePThreeCFluxVariables<TypeTag>);

//! the BoundaryVariables property
// SET_TYPE_PROP(BoxHolleThreePThreeC, BoundaryVariables, HolleThreePThreeCBoundaryVariables<TypeTag>);

//! the upwind factor for the mobility.
SET_SCALAR_PROP(BoxHolleThreePThreeC, MobilityUpwindAlpha, 1.0);

//! The indices required by the isothermal holle3p3c model
SET_PROP(BoxHolleThreePThreeC,
         HolleThreePThreeCIndices)
{ private:
    enum { Formulation = GET_PROP_VALUE(TypeTag, PTAG(Formulation)) };
public:
    typedef HolleThreePThreeCIndices<TypeTag, Formulation, 0> type;
};

// enable jacobian matrix recycling by default
SET_BOOL_PROP(BoxHolleThreePThreeC, EnableJacobianRecycling, true);
// enable partial reassembling by default
SET_BOOL_PROP(BoxHolleThreePThreeC, EnablePartialReassemble, true);
// enable time-step ramp up by default
SET_BOOL_PROP(BoxHolleThreePThreeC, EnableTimeStepRampUp, false);

//
}

}

#endif
