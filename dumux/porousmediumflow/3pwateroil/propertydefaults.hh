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
 * \ingroup ThreePWaterOilModel
 */
/*!
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        3p2cni fully implicit model.
 */
#ifndef DUMUX_3P2CNI_PROPERTY_DEFAULTS_HH
#define DUMUX_3P2CNI_PROPERTY_DEFAULTS_HH

#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>

#include "indices.hh"
#include "model.hh"
#include "fluxvariables.hh"
#include "volumevariables.hh"
#include "properties.hh"
#include "newtoncontroller.hh"
#include "localresidual.hh"

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
SET_PROP(ThreePWaterOilNI, NumComponents)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

 public:
    static const int value = FluidSystem::numComponents;

    static_assert(value == 2,
                  "Only fluid systems with 2 components are supported by the 3p2cni model!");
};

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(ThreePWaterOilNI, NumPhases)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

 public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 3,
                  "Only fluid systems with 3 phases are supported by the 3p2cni model!");
};

//! Use the 3p2cni specific newton controller for the 3p2cni model
SET_TYPE_PROP(ThreePWaterOilNI, NewtonController, ThreePWaterOilNewtonController<TypeTag>);

//! define the base flux variables to realize Darcy flow
SET_TYPE_PROP(ThreePWaterOilNI, BaseFluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//! the upwind factor for the mobility.
SET_SCALAR_PROP(ThreePWaterOilNI, ImplicitMassUpwindWeight, 1.0);

//! set default mobility upwind weight to 1.0, i.e. fully upwind
SET_SCALAR_PROP(ThreePWaterOilNI, ImplicitMobilityUpwindWeight, 1.0);

//! Determines whether a constraint solver should be used explicitly
SET_BOOL_PROP(ThreePWaterOilNI, UseSimpleModel, true);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(ThreePWaterOilNI, SpatialParams, ImplicitSpatialParams<TypeTag>);

//! Use the model after Millington (1961) for the effective diffusivity
SET_TYPE_PROP(ThreePWaterOilNI, EffectiveDiffusivityModel,
             DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// disable velocity output by default

// enable gravity by default
SET_BOOL_PROP(ThreePWaterOilNI, ProblemEnableGravity, true);

SET_BOOL_PROP(ThreePWaterOilNI, UseMoles, true); //!< Define that mole fractions are used in the balance equations per default



//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(ThreePWaterOilNI, SpatialParamsForchCoeff, 0.55);


//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(ThreePWaterOilNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

//! temperature is already written by the isothermal model
SET_BOOL_PROP(ThreePWaterOilNI, NiOutputLevel, 0);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(ThreePWaterOilNI, IsothermalModel, ThreePWaterOilModel<TypeTag>);

// set isothermal FluxVariables
SET_TYPE_PROP(ThreePWaterOilNI, IsothermalFluxVariables, ThreePWaterOilFluxVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(ThreePWaterOilNI, IsothermalVolumeVariables, ThreePWaterOilVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(ThreePWaterOilNI, IsothermalLocalResidual, ThreePWaterOilLocalResidual<TypeTag>);

//set isothermal Indices
SET_PROP(ThreePWaterOilNI, IsothermalIndices)
{

public:
    typedef ThreePWaterOilIndices<TypeTag, /*PVOffset=*/0> type;
};

//set isothermal NumEq
SET_INT_PROP(ThreePWaterOilNI, IsothermalNumEq, 2);

}

}

#endif
