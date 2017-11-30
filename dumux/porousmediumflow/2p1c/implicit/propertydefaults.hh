// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \ingroup TwoPOneCModel
 */
/*!
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        2p1cni fully implicit model.
 *
 *Important note: The 2p1c model requires the use of the non-isothermal extension found in dumux/implicit/nonisothermal
 */
#ifndef DUMUX_2P1C_PROPERTY_DEFAULTS_HH
#define DUMUX_2P1C_PROPERTY_DEFAULTS_HH

#include "model.hh"
#include "indices.hh"
#include "volumevariables.hh"
#include "properties.hh"
#include "newtoncontroller.hh"
#include "localresidual.hh"

#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/fluidsystems/2pliquidvapor.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/spatialparams/implicit.hh>

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
 * assert to make sure it is 1.
 */
SET_PROP(TwoPOneCNI, NumComponents)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

 public:
    static const int value = FluidSystem::numComponents;

    static_assert(value == 1,
                  "Only fluid systems with 1 component are supported by the 2p1cni model!");
};

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(TwoPOneCNI, NumPhases)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

 public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "Only fluid systems with 2 phases are supported by the 2p1cni model!");
};

//! Use the 2p2c specific newton controller for dealing with phase switches
SET_TYPE_PROP(TwoPOneCNI, NewtonController, TwoPOneCNINewtonController<TypeTag>);

//! the FluxVariables property
// SET_TYPE_PROP(TwoPOneCNI, FluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//! the upwind factor for the mobility.
SET_SCALAR_PROP(TwoPOneCNI, ImplicitMassUpwindWeight, 1.0);

//! set default mobility upwind weight to 1.0, i.e. fully upwind
SET_SCALAR_PROP(TwoPOneCNI, ImplicitMobilityUpwindWeight, 1.0);

//! The fluid system to use by default
SET_PROP(TwoPOneCNI, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::TwoPLiquidVaporFluidsystem<Scalar, NullComponent<Scalar> > type;
};

//! Determines whether Blocking ofspurious flow is used
SET_BOOL_PROP(TwoPOneCNI, UseBlockingOfSpuriousFlow, false);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(TwoPOneCNI, SpatialParams, ImplicitSpatialParams<TypeTag>);

// disable velocity output by default

// enable gravity by default
SET_BOOL_PROP(TwoPOneCNI, ProblemEnableGravity, true);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(BoxModel, SpatialParamsForchCoeff, 0.55);


//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(NonIsothermal, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar> type;
};

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(TwoPOneCNI, IsothermalModel, TwoPOneCModel<TypeTag>);

// set isothermal FluxVariables
SET_TYPE_PROP(TwoPOneCNI, IsothermalFluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(TwoPOneCNI, IsothermalVolumeVariables, TwoPOneCVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(TwoPOneCNI, IsothermalLocalResidual, TwoPOneCLocalResidual<TypeTag>);

//set isothermal Indices
SET_PROP(TwoPOneCNI, IsothermalIndices)
{

public:
    typedef TwoPOneCIndices<TypeTag, 0> type;
};

//set isothermal NumEq
SET_INT_PROP(TwoPOneCNI, IsothermalNumEq, 1);

}

}

#endif
