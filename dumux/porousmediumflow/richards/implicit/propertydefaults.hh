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
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup RichardsModel
 * \file
 *
 * \brief Contains the default definitions for the properties required
 *        by the Richards fully implicit model.
 */
#ifndef DUMUX_RICHARDS_PROPERTY_DEFAULTS_HH
#define DUMUX_RICHARDS_PROPERTY_DEFAULTS_HH

#include "model.hh"
#include "indices.hh"
#include "volumevariables.hh"
#include "properties.hh"
#include "newtoncontroller.hh"
#include "localresidual.hh"
#include "primaryvariableswitch.hh"

#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidstates/immiscible.hh>

namespace Dumux
{
// \{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties values
//////////////////////////////////////////////////////////////////
//! Number of equations required by the model
SET_INT_PROP(Richards, NumEq, 1);

//! Number of fluid phases considered
SET_INT_PROP(Richards, NumPhases, 2);

//! Number of components considered (only water)
SET_INT_PROP(Richards, NumComponents, 1);

//! The local residual operator
SET_TYPE_PROP(Richards, LocalResidual, RichardsLocalResidual<TypeTag>);

//! The global model used
SET_TYPE_PROP(Richards, Model, RichardsModel<TypeTag>);

//! The class for the volume averaged quantities
SET_TYPE_PROP(Richards, VolumeVariables, RichardsVolumeVariables<TypeTag>);

//! Smarter newton controller
SET_TYPE_PROP(Richards, NewtonController, RichardsNewtonController<TypeTag>);

//! Enable advection
SET_BOOL_PROP(Richards, EnableAdvection, true);

//! The default richards model computes no diffusion in the air phase
//! Turning this on leads to the extended Richards equation (see e.g. Vanderborght et al. 2017)
SET_BOOL_PROP(Richards, EnableWaterDiffusionInAir, false);

//! The default richards model has no molecular diffusion
SET_BOOL_PROP(Richards, EnableMolecularDiffusion, GET_PROP_VALUE(TypeTag, EnableWaterDiffusionInAir));

//! Use the model after Millington (1961) for the effective diffusivity
SET_TYPE_PROP(Richards, EffectiveDiffusivityModel,
              DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! The default is not to use the kelvin equation for the water vapor pressure (dependency on pc)
SET_BOOL_PROP(Richards, UseKelvinEquation, false);

//! Isothermal model by default
SET_BOOL_PROP(Richards, EnableEnergyBalance, false);

//! The class with all index definitions for the model
SET_TYPE_PROP(Richards, Indices, RichardsIndices);

//! The class with all index definitions for the model
SET_TYPE_PROP(Richards, PrimaryVariables, SwitchablePrimaryVariables<TypeTag, int>);

//! The primary variable switch for the richards model
SET_TYPE_PROP(Richards, PrimaryVariableSwitch, ExtendedRichardsPrimaryVariableSwitch<TypeTag>);

//! The primary variable switch for the richards model
SET_BOOL_PROP(Richards, ProblemUsePrimaryVariableSwitch, false);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(Richards, SpatialParams, ImplicitSpatialParams<TypeTag>);

/*!
 * \brief Set type of the parameter objects for the material law
 *
 * By default this is just retrieved from the material law.
 */
SET_TYPE_PROP(Richards, MaterialLawParams, typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

/*!
 *\brief The fluid system used by the model.
 *
 * By default this uses the H2O-Air fluid system with Simple H2O (constant density and viscosity).
 */
SET_PROP(Richards, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::H2OAir<Scalar, SimpleH2O<Scalar>, false>;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(Richards, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = ImmiscibleFluidState<Scalar, FluidSystem>;
};

// enable gravity by default
SET_BOOL_PROP(Richards, ProblemEnableGravity, true);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(Richards, SpatialParamsForchCoeff, 0.55);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(RichardsNI, ThermalConductivityModel)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
public:
    using type = ThermalConductivitySomerton<Scalar, Indices>;
};

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(RichardsNI, IsothermalModel, RichardsModel<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(RichardsNI, IsothermalVolumeVariables, RichardsVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(RichardsNI, IsothermalLocalResidual, ImmiscibleLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(RichardsNI, IsothermalIndices, RichardsIndices);

//set isothermal NumEq
SET_INT_PROP(RichardsNI, IsothermalNumEq, 1);

// \}
}

} // end namespace

#endif
