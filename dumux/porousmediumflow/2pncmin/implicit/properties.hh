// -**- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \ingroup TwoPNCMinModel
 *
 * \file
 *
 * \brief Defines the properties required for the two-phase n-component mineralization
 *        fully implicit model.
 */
#ifndef DUMUX_2PNCMIN_PROPERTIES_HH
#define DUMUX_2PNCMIN_PROPERTIES_HH

#include <dumux/porousmediumflow/2pnc/implicit/properties.hh>

#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
NEW_TYPE_TAG(TwoPNCMin, INHERITS_FROM(TwoPNC));
NEW_TYPE_TAG(TwoPNCMinNI, INHERITS_FROM(TwoPNCMin, NonIsothermal));

//////////////////////////////////////////////////////////////////
// Property tags for the isothermal 2pncmin model
//////////////////////////////////////////////////////////////////

SET_TYPE_PROP(TwoPNCMin, PrimaryVariables, SwitchablePrimaryVariables<TypeTag, int>);          //! The primary variables vector for the 2pncmin model
SET_TYPE_PROP(TwoPNCMin, PrimaryVariableSwitch, TwoPNCPrimaryVariableSwitch<TypeTag>);         //! Using the 2pnc primary variable switch for the 2pncmin model
SET_TYPE_PROP(TwoPNCMin, VolumeVariables, TwoPNCMinVolumeVariables<TypeTag>);                  //! the VolumeVariables property
SET_TYPE_PROP(TwoPNCMin, Indices, TwoPNCIndices <TypeTag, /*PVOffset=*/0>);                    //! Using the 2pnc indices required by the isothermal 2pncmin model
SET_TYPE_PROP(TwoPNCMin, SpatialParams, ImplicitSpatialParams<TypeTag>);                       //! Use the ImplicitSpatialParams by default
SET_TYPE_PROP(TwoPNCMin, VtkOutputFields, TwoPNCMinVtkOutputFields<TypeTag>);                  //! Set the vtk output fields specific to the TwoPNCMin model
SET_TYPE_PROP(TwoPNCMin, LocalResidual, CompositionalLocalResidual<TypeTag>);                  //! Use the compositional local residual

SET_INT_PROP(TwoPNCMin, NumComponents, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents);    //! Use the number of components of the fluid system
SET_INT_PROP(TwoPNCMin, ReplaceCompEqIdx, GET_PROP_TYPE(TypeTag, FluidSystem)::numComponents); //! Per default, no component mass balance is replaced
SET_INT_PROP(TwoPNCMin, Formulation, TwoPNCFormulation::pwsn);                                 //! Using the default 2pnc formulation: pw-Sn, overwrite if necessary

SET_BOOL_PROP(TwoPNCMin, SetMoleFractionsForWettingPhase, true);  //! Set the primary variables mole fractions for the wetting or non-wetting phase
SET_BOOL_PROP(TwoPNCMin, EnableAdvection, true);                  //! Enable advection
SET_BOOL_PROP(TwoPNCMin, EnableMolecularDiffusion, true);         //! Enable molecular diffusion
SET_BOOL_PROP(TwoPNCMin, EnableEnergyBalance, false);             //! This is the isothermal variant of the model
SET_BOOL_PROP(TwoPNCMin, UseMoles, true);                         //! Use mole fractions in the balance equations by default

//! Use the model after Millington (1961) for the effective diffusivity
SET_TYPE_PROP(TwoPNCMin, EffectiveDiffusivityModel, DiffusivityMillingtonQuirk<typename GET_PROP_TYPE(TypeTag, Scalar)>);

//! The major components belonging to the existing phases, e.g. 2 for water and air being the major components in a liquid-gas-phase system
SET_PROP(TwoPNCMin, NumMajorComponents)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2, "The model is restricted to two phases, thus number of major components must also be two.");
};

//! We use the number of phases of the fluid system. Make sure it is 2!
SET_PROP(TwoPNCMin, NumPhases)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2, "Only fluid systems with 2 fluid phases are supported by the 2p-nc-min model!");
};

//! This model uses the compositional fluid state
SET_PROP(TwoPNCMin, FluidState)
{
private:
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
public:
    using type = CompositionalFluidState<Scalar, FluidSystem>;
};

//! Set the property for the material parameters by extracting it from the material law.
SET_PROP(TwoPNCMin, MaterialLawParams)
{
private:
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw));

public:
    using type = typename MaterialLaw::Params;
};

//! Set the property for the number of solid phases, excluding the non-reactive matrix.
SET_PROP(TwoPNCMin, NumSPhases)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));

public:
    static const int value = FluidSystem::numSPhases;
};

//! Set the property for the number of equations. For each component and each
//precipitated mineral/solid phase one equation has to be solved.

SET_PROP(TwoPNCMin, NumEq)
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem));

public:
    static const int value = FluidSystem::numComponents + FluidSystem::numSPhases;
};
/////////////////////////////////////////////////
// Properties for the non-isothermal 2pncmin model
/////////////////////////////////////////////////
SET_TYPE_PROP(TwoPNCMinNI, IsothermalVolumeVariables, TwoPNCMinVolumeVariables<TypeTag>);     //! set isothermal VolumeVariables
SET_TYPE_PROP(TwoPNCMinNI, IsothermalLocalResidual, CompositionalLocalResidual<TypeTag>);  //! set isothermal LocalResidual
SET_TYPE_PROP(TwoPNCMinNI, IsothermalIndices, TwoPNCMinIndices<TypeTag, /*PVOffset=*/0>);     //! set isothermal Indices
SET_TYPE_PROP(TwoPNCMinNI, IsothermalVtkOutputFields, TwoPNCMinVtkOutputFields<TypeTag>);     //! set isothermal output fields

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(TwoPNCMinNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

//! set isothermal NumEq
SET_PROP(TwoPNCMinNI, IsothermalNumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;
};

}
}

#endif
