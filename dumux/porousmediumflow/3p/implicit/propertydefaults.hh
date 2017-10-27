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

#include "indices.hh"

#include "model.hh"
#include "indices.hh"
#include "volumevariables.hh"
#include "properties.hh"
#include "vtkoutputfields.hh"

#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidstates/immiscible.hh>
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

SET_PROP(ThreeP, NumComponents)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

 public:
    static const int value = FluidSystem::numComponents;
    static_assert(value == 3,
                  "Only fluid systems with 3 components are supported by the 3p model!");
};

SET_INT_PROP(ThreeP, NumEq, 3); //!< set the number of equations to 2

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(ThreeP, MaterialLawParams, typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! The local residual function of the conservation equations
SET_TYPE_PROP(ThreeP, LocalResidual, ImmiscibleLocalResidual<TypeTag>);

//! Enable advection
SET_BOOL_PROP(ThreeP, EnableAdvection, true);

//! disable molecular diffusion for the 3p model
SET_BOOL_PROP(ThreeP, EnableMolecularDiffusion, false);

//! Isothermal model by default
SET_BOOL_PROP(ThreeP, EnableEnergyBalance, false);

//! the Model property
SET_TYPE_PROP(ThreeP, Model, ThreePModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(ThreeP, VolumeVariables, ThreePVolumeVariables<TypeTag>);

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

// enable gravity by default
SET_BOOL_PROP(ThreeP, ProblemEnableGravity, true);

//! Set the vtk output fields specific to the ThreeP model
SET_TYPE_PROP(ThreeP, VtkOutputFields, ThreePVtkOutputFields<TypeTag>);

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

//set isothermal VolumeVariables
SET_TYPE_PROP(ThreePNI, IsothermalVolumeVariables, ThreePVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(ThreePNI, IsothermalLocalResidual, ImmiscibleLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(ThreePNI, IsothermalIndices, ThreePIndices<TypeTag,/*PVOffset=*/0>);

//set isothermal NumEq
SET_INT_PROP(ThreePNI, IsothermalNumEq, 3);

}

}

#endif
