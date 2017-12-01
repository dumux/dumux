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
 * \ingroup ThreePThreeCModel
 */
/*!
 * \file
 *
 * \brief Defines the properties required for the three-phase three-component
 *        fully implicit model.
 */
#ifndef DUMUX_3P3C_PROPERTIES_HH
#define DUMUX_3P3C_PROPERTIES_HH

#include <dumux/common/properties.hh>

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/properties.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"
#include "primaryvariableswitch.hh"
#include "localresidual.hh"

#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>

namespace Dumux
{

namespace Properties
{
//! The type tags for the implicit three-phase three-component problems
NEW_TYPE_TAG(ThreePThreeC, INHERITS_FROM(PorousMediumFlow));

//! The type tags for the corresponding non-isothermal problems
NEW_TYPE_TAG(ThreePThreeCNI, INHERITS_FROM(ThreePThreeC, NonIsothermal));

//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of components.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 3.
 */
SET_PROP(ThreePThreeC, NumComponents)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

 public:
    static const int value = FluidSystem::numComponents;

    static_assert(value == 3,
                  "Only fluid systems with 3 components are supported by the 3p3c model!");
};

/*!
 * \brief Set the property for the number of fluid phases.
 *
 * We just forward the number from the fluid system and use an static
 * assert to make sure it is 3.
 */
SET_PROP(ThreePThreeC, NumPhases)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

 public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 3,
                  "Only fluid systems with 3 phases are supported by the 3p3c model!");
};

//! Set as default that no component mass balance is replaced by the total mass balance
SET_INT_PROP(ThreePThreeC, ReplaceCompEqIdx, 100);
/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(ThreePThreeC, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef CompositionalFluidState<Scalar, FluidSystem> type;
};

SET_INT_PROP(ThreePThreeC, NumEq, 3); //!< set the number of equations to 2

//! The local residual function of the conservation equations
SET_TYPE_PROP(ThreePThreeC, LocalResidual, ThreePThreeCLocalResidual<TypeTag>);

//! Enable advection
SET_BOOL_PROP(ThreePThreeC, EnableAdvection, true);

//! Enable molecular diffusion
SET_BOOL_PROP(ThreePThreeC, EnableMolecularDiffusion, true);

//! Isothermal model by default
SET_BOOL_PROP(ThreePThreeC, EnableEnergyBalance, false);

//! The primary variable switch for the 3p3c model
SET_TYPE_PROP(ThreePThreeC, PrimaryVariableSwitch, ThreePThreeCPrimaryVariableSwitch<TypeTag>);

//! The primary variables vector for the 2p2c model
SET_TYPE_PROP(ThreePThreeC, PrimaryVariables, SwitchablePrimaryVariables<TypeTag, int>);

//! the VolumeVariables property
SET_TYPE_PROP(ThreePThreeC, VolumeVariables, ThreePThreeCVolumeVariables<TypeTag>);

//! Determines whether a constraint solver should be used explicitly
SET_BOOL_PROP(ThreePThreeC, UseConstraintSolver, false);

//! The indices required by the isothermal 3p3c model
SET_TYPE_PROP(ThreePThreeC, Indices, ThreePThreeCIndices<TypeTag, /*PVOffset=*/0>);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(ThreePThreeC, SpatialParams, ImplicitSpatialParams<TypeTag>);

//! The model after Millington (1961) is used for the effective diffusivity
SET_PROP(ThreePThreeC, EffectiveDiffusivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef DiffusivityMillingtonQuirk<Scalar> type;
};

//! Set the vtk output fields specific to the ThreeP model
SET_TYPE_PROP(ThreePThreeC, VtkOutputFields, ThreePThreeCVtkOutputFields<TypeTag>);

//! Use mole fractions in the balance equations by default
SET_BOOL_PROP(ThreePThreeC, UseMoles, true);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(ThreePThreeCNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

//set isothermal VolumeVariables
SET_TYPE_PROP(ThreePThreeCNI, IsothermalVolumeVariables, ThreePThreeCVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(ThreePThreeCNI, IsothermalLocalResidual, ThreePThreeCLocalResidual<TypeTag>);

//set isothermal Indices
SET_PROP(ThreePThreeCNI, IsothermalIndices)
{

public:
    typedef ThreePThreeCIndices<TypeTag, /*PVOffset=*/0> type;
};

//set isothermal NumEq
SET_INT_PROP(ThreePThreeCNI, IsothermalNumEq, 3);

//! Set the vtk output fields specific to the ThreeP model
SET_TYPE_PROP(ThreePThreeCNI, IsothermalVtkOutputFields, ThreePThreeCVtkOutputFields<TypeTag>);

}

}

#endif
