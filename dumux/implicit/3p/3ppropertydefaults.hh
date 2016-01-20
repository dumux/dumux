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
 * \ingroup ThreePModel
 */
/*!
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        three-phase fully implicit model.
 */
#ifndef DUMUX_3P_PROPERTY_DEFAULTS_HH
#define DUMUX_3P_PROPERTY_DEFAULTS_HH

#include "3pindices.hh"

#include "3pmodel.hh"
#include "3pindices.hh"
#include "3pvolumevariables.hh"
#include "3pproperties.hh"
#include "3plocalresidual.hh"

#include <dumux/implicit/nonisothermal/nipropertydefaults.hh>
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>


namespace Dumux
{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(ThreeP, NumPhases)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

 public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 3,
                  "Only fluid systems with 3 phases are supported by the 3p model!");
};

SET_INT_PROP(ThreeP, NumEq, 3); //!< set the number of equations to 2

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(ThreeP, MaterialLawParams, typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! The local residual function of the conservation equations
SET_TYPE_PROP(ThreeP, LocalResidual, ThreePLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(ThreeP, Model, ThreePModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(ThreeP, VolumeVariables, ThreePVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(ThreeP, FluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//! the upwind factor for the mobility.
SET_SCALAR_PROP(ThreeP, ImplicitMassUpwindWeight, 1.0);

//! set default mobility upwind weight to 1.0, i.e. fully upwind
SET_SCALAR_PROP(ThreeP, ImplicitMobilityUpwindWeight, 1.0);

//! The indices required by the isothermal 3p model
SET_TYPE_PROP(ThreeP, Indices, ThreePIndices<TypeTag,/*PVOffset=*/0>);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(ThreeP, SpatialParams, ImplicitSpatialParams<TypeTag>);

SET_PROP(ThreeP, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};

// disable velocity output by default
SET_BOOL_PROP(ThreeP, VtkAddVelocity, false);

// enable gravity by default
SET_BOOL_PROP(ThreeP, ProblemEnableGravity, true);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(ThreePNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

//! temperature is already written by the isothermal model
SET_BOOL_PROP(ThreePNI, NiOutputLevel, 0);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(ThreePNI, IsothermalModel, ThreePModel<TypeTag>);

// set isothermal FluxVariables
SET_TYPE_PROP(ThreePNI, IsothermalFluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(ThreePNI, IsothermalVolumeVariables, ThreePVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(ThreePNI, IsothermalLocalResidual, ThreePLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(ThreePNI, IsothermalIndices, ThreePIndices<TypeTag,/*PVOffset=*/0>);

//set isothermal NumEq
SET_INT_PROP(ThreePNI, IsothermalNumEq, 3);

}

}

#endif
