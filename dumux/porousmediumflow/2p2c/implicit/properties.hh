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
 * \ingroup TwoPTwoCModel
 *
 * \file
 *
 * \brief Defines the properties required for the two-phase two-component
 *        fully implicit model.
 */
#ifndef DUMUX_2P2C_PROPERTIES_HH
#define DUMUX_2P2C_PROPERTIES_HH

// property forward declarations
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/model.hh>

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>

#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/porousmediumflow/compositional/privarswitchnewtoncontroller.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "primaryvariableswitch.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
NEW_TYPE_TAG(TwoPTwoC, INHERITS_FROM(PorousMediumFlow));
NEW_TYPE_TAG(TwoPTwoCNI, INHERITS_FROM(TwoPTwoC, NonIsothermal));

//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

//! Set the number of equations to 2
SET_INT_PROP(TwoPTwoC, NumEq, 2);

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system and use a static
 * assert to make sure it is 2.
 */
SET_PROP(TwoPTwoC, NumComponents)
{
    static constexpr int value = 2;
    static_assert(GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents == value,
                  "Only fluid systems with 2 components are supported by the 2p-2c model!");
};

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use a static
 * assert to make sure it is 2.
 */
SET_PROP(TwoPTwoC, NumPhases)
{
    static constexpr int value = 2;
    static_assert(GET_PROP_TYPE(TypeTag, FluidSystem)::numPhases == value,
                  "Only fluid systems with 2 phases are supported by the 2p-2c model!");
};

//! Set the vtk output fields specific to the TwoPTwoC model
SET_TYPE_PROP(TwoPTwoC, VtkOutputFields, TwoPTwoCVtkOutputFields<TypeTag>);

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(TwoPTwoC, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! Set the default formulation to pw-sn
SET_INT_PROP(TwoPTwoC, Formulation, TwoPTwoCFormulation::pwsn);

//! Set as default that no component mass balance is replaced by the total mass balance
SET_INT_PROP(TwoPTwoC, ReplaceCompEqIdx, GET_PROP_VALUE(TypeTag, NumComponents));

//! Use the compositional local residual operator
SET_TYPE_PROP(TwoPTwoC, LocalResidual, CompositionalLocalResidual<TypeTag>);

//! Enable advection
SET_BOOL_PROP(TwoPTwoC, EnableAdvection, true);

//! Enable molecular diffusion
SET_BOOL_PROP(TwoPTwoC, EnableMolecularDiffusion, true);

//! Isothermal model by default
SET_BOOL_PROP(TwoPTwoC, EnableEnergyBalance, false);

//! The primary variable switch for the 2p2c model
SET_TYPE_PROP(TwoPTwoC, PrimaryVariableSwitch, TwoPTwoCPrimaryVariableSwitch<TypeTag>);

//! The primary variables vector for the 2p2c model
SET_TYPE_PROP(TwoPTwoC, PrimaryVariables, SwitchablePrimaryVariables<TypeTag, int>);

//! Use the 2p2c VolumeVariables
SET_TYPE_PROP(TwoPTwoC, VolumeVariables, TwoPTwoCVolumeVariables<TypeTag>);

//! Set the indices required by the isothermal 2p2c
SET_TYPE_PROP(TwoPTwoC, Indices, TwoPTwoCIndices<typename GET_PROP_TYPE(TypeTag, FluidSystem), /*PVOffset=*/0>);

//! Use the ImplicitSpatialParams by default
SET_TYPE_PROP(TwoPTwoC, SpatialParams, ImplicitSpatialParams<TypeTag>);

//! Use the model after Millington (1961) for the effective diffusivity
SET_TYPE_PROP(TwoPTwoC, EffectiveDiffusivityModel,
             DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! Use mole fractions in the balance equations by default
SET_BOOL_PROP(TwoPTwoC, UseMoles, true);

//! Determines whether the constraint solver is used
SET_BOOL_PROP(TwoPTwoC, UseConstraintSolver, true);

//! Determines whether the Kelvin equation is used to adapt the saturation vapor pressure
SET_BOOL_PROP(TwoPTwoC, UseKelvinEquation, false);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(TwoPTwoCNI, ThermalConductivityModel)
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

//set isothermal Indices
SET_TYPE_PROP(TwoPTwoCNI, IsothermalIndices, TwoPTwoCIndices<typename GET_PROP_TYPE(TypeTag, FluidSystem), /*PVOffset=*/0>);

//set isothermal output fields
SET_TYPE_PROP(TwoPTwoCNI, IsothermalVtkOutputFields, TwoPTwoCVtkOutputFields<TypeTag>);

//set isothermal NumEq
SET_INT_PROP(TwoPTwoCNI, IsothermalNumEq, 2);

} // end namespace Properties
} // end namespace Dumux

#endif
