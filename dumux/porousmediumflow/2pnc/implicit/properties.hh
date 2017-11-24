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
 *
 * \file
 *
 * \brief Defines the properties required for the two-phase n-component
 *        fully implicit model.
 */
#ifndef DUMUX_2PNC_PROPERTIES_HH
#define DUMUX_2PNC_PROPERTIES_HH

#include <dumux/common/basicproperties.hh>
#include <dumux/linear/linearsolverproperties.hh>

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/compositional/switchableprimaryvariables.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/properties.hh>

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
NEW_TYPE_TAG(TwoPNC, INHERITS_FROM(PorousMediumFlow, NumericModel, LinearSolverTypeTag));
NEW_TYPE_TAG(TwoPNCNI, INHERITS_FROM(TwoPNC, NonIsothermal));

//////////////////////////////////////////////////////////////////
// Properties for the isothermal 2pnc model
//////////////////////////////////////////////////////////////////
SET_TYPE_PROP(TwoPNC, PrimaryVariables, SwitchablePrimaryVariables<TypeTag, int>);          //! The primary variables vector for the 2pnc model
SET_TYPE_PROP(TwoPNC, PrimaryVariableSwitch, TwoPNCPrimaryVariableSwitch<TypeTag>);         //! The primary variable switch for the 2pnc model
SET_TYPE_PROP(TwoPNC, VolumeVariables, TwoPNCVolumeVariables<TypeTag>);                     //! the VolumeVariables property
SET_TYPE_PROP(TwoPNC, Indices, TwoPNCIndices <TypeTag, /*PVOffset=*/0>);                    //! The indices required by the isothermal 2pnc model
SET_TYPE_PROP(TwoPNC, SpatialParams, ImplicitSpatialParams<TypeTag>);                       //! Use the ImplicitSpatialParams by default
SET_TYPE_PROP(TwoPNC, VtkOutputFields, TwoPNCVtkOutputFields<TypeTag>);                     //! Set the vtk output fields specific to the TwoPNC model
SET_TYPE_PROP(TwoPNC, LocalResidual, CompositionalLocalResidual<TypeTag>);                  //! Use the compositional local residual

SET_INT_PROP(TwoPNC, NumComponents, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents);    //! Use the number of components of the fluid system
SET_INT_PROP(TwoPNC, ReplaceCompEqIdx, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents); //! Per default, no component mass balance is replaced
SET_INT_PROP(TwoPNC, NumEq, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents);            //! We solve one equation per component
SET_INT_PROP(TwoPNC, Formulation, TwoPNCFormulation::pwsn);                                 //! Default formulation is pw-Sn, overwrite if necessary

SET_BOOL_PROP(TwoPNC, SetMoleFractionsForWettingPhase, true);  //! Set the primary variables mole fractions for the wetting or non-wetting phase
SET_BOOL_PROP(TwoPNC, EnableAdvection, true);                  //! Enable advection
SET_BOOL_PROP(TwoPNC, EnableMolecularDiffusion, true);         //! Enable molecular diffusion
SET_BOOL_PROP(TwoPNC, EnableEnergyBalance, false);             //! This is the isothermal variant of the model
SET_BOOL_PROP(TwoPNC, UseMoles, true);                         //! Use mole fractions in the balance equations by default


//! Use the model after Millington (1961) for the effective diffusivity
SET_TYPE_PROP(TwoPNC, EffectiveDiffusivityModel, DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! The major components belonging to the existing phases, e.g. 2 for water and air being the major components in a liquid-gas-phase system
SET_PROP(TwoPNC, NumMajorComponents)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2, "The model is restricted to two phases, thus number of major components must also be two.");
};

//! We use the number of phases of the fluid system. Make sure it is 2!
SET_PROP(TwoPNC, NumPhases)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2, "Only fluid systems with 2 fluid phases are supported by the 2p-nc model!");
};

//! This model uses the compositional fluid state
SET_PROP(TwoPNC, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! Set the property for the material parameters by extracting it from the material law.
SET_PROP(TwoPNC, MaterialLawParams)
{
private:
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw));

public:
    using type = typename MaterialLaw::Params;
};

/////////////////////////////////////////////////
// Properties for the non-isothermal 2pnc model
/////////////////////////////////////////////////
SET_TYPE_PROP(TwoPNCNI, IsothermalVolumeVariables, TwoPNCVolumeVariables<TypeTag>);    //! set isothermal VolumeVariables
SET_TYPE_PROP(TwoPNCNI, IsothermalLocalResidual, CompositionalLocalResidual<TypeTag>); //! set isothermal LocalResidual
SET_TYPE_PROP(TwoPNCNI, IsothermalIndices, TwoPNCIndices<TypeTag, /*PVOffset=*/0>);    //! set isothermal Indices

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(TwoPNCNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

//! set isothermal NumEq
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
