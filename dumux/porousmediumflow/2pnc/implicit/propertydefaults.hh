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
 * \ingroup TwoPNCModel
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        two-phase n-component fully implicit model.
 */
#ifndef DUMUX_2PNC_PROPERTY_DEFAULTS_HH
#define DUMUX_2PNC_PROPERTY_DEFAULTS_HH

#include "indices.hh"
#include "model.hh"
#include "volumevariables.hh"
#include "properties.hh"
#include "newtoncontroller.hh"
#include "primaryvariableswitch.hh"

#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>

namespace Dumux
{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system
 *
 */
SET_PROP(TwoPNC, NumComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;

};
//! Set as default that no component mass balance is replaced by the total mass balance
SET_PROP(TwoPNC, ReplaceCompEqIdx)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;
};
//! The major components belonging to the existing phases, e.g. 2 for water and air being the major components in a liquid-gas-phase system
SET_PROP(TwoPNC, NumMajorComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "The model is restricted to two phases, thus number of major components must also be two.");
};

//! Set the primary variables mole fractions for the wetting or non-wetting phase
SET_BOOL_PROP(TwoPNC, SetMoleFractionsForWettingPhase, true);

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 2.
 */
SET_PROP(TwoPNC, NumPhases)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "Only fluid systems with 2 fluid phases are supported by the 2p-nc model!");
};
/*!
 * \brief Set the property for the number of equations: For each existing component one equation has to be solved.
 */
SET_PROP(TwoPNC, NumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(TwoPNC, FluidState)
{
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef CompositionalFluidState<Scalar, FluidSystem> type;
};

//! Set the default formulation to pw-Sn: This can be over written in the problem.
SET_INT_PROP(TwoPNC, Formulation, TwoPNCFormulation::pwsn);

//! Set the property for the material parameters by extracting it from the material law.
SET_PROP(TwoPNC, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

//! Use the 2pnc local residual
SET_TYPE_PROP(TwoPNC, LocalResidual, CompositionalLocalResidual<TypeTag>);

//! Enable advection
SET_BOOL_PROP(TwoPNC, EnableAdvection, true);

//! Enable molecular diffusion
SET_BOOL_PROP(TwoPNC, EnableMolecularDiffusion, true);

//! Isothermal model by default
SET_BOOL_PROP(TwoPNC, EnableEnergyBalance, false);

//! Use the 2pnc newton controller
SET_TYPE_PROP(TwoPNC, NewtonController, TwoPNCNewtonController<TypeTag>);

//! the Model property
SET_TYPE_PROP(TwoPNC, Model, TwoPNCModel<TypeTag>);

//! The primary variable switch for the 2pnc model
SET_TYPE_PROP(TwoPNC, PrimaryVariableSwitch, TwoPNCPrimaryVariableSwitch<TypeTag>);

//! The primary variables vector for the 2pnc model
SET_TYPE_PROP(TwoPNC, PrimaryVariables, SwitchablePrimaryVariables<TypeTag, int>);

//! the VolumeVariables property
SET_TYPE_PROP(TwoPNC, VolumeVariables, TwoPNCVolumeVariables<TypeTag>);

//! The indices required by the isothermal 2pnc model
SET_TYPE_PROP(TwoPNC, Indices, TwoPNCIndices <TypeTag, /*PVOffset=*/0>);

//! Use the ImplicitSpatialParams by default
SET_TYPE_PROP(TwoPNC, SpatialParams, ImplicitSpatialParams<TypeTag>);

//! Use the model after Millington (1961) for the effective diffusivity
SET_TYPE_PROP(TwoPNC, EffectiveDiffusivityModel, DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! Enable gravity by default
SET_BOOL_PROP(TwoPNC, ProblemEnableGravity, true);

//! Use mole fractions in the balance equations by default
SET_BOOL_PROP(TwoPNC, UseMoles, true);

//! Disable velocity output by default

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(TwoPNCNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

//! temperature is already written by the isothermal model
SET_BOOL_PROP(TwoPNCNI, NiOutputLevel, 0);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(TwoPNCNI, IsothermalModel, TwoPNCModel<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(TwoPNCNI, IsothermalVolumeVariables, TwoPNCVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(TwoPNCNI, IsothermalLocalResidual, CompositionalLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(TwoPNCNI, IsothermalIndices, TwoPNCIndices<TypeTag, /*PVOffset=*/0>);

//set isothermal NumEq
SET_PROP(TwoPNCNI, IsothermalNumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;
};


}

}

#endif
