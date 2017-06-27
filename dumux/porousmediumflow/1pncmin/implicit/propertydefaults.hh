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
 * \ingroup OnePNCMinModel
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        two-phase n-component mineralization fully implicit model.
 */
#ifndef DUMUX_1PNCMIN_PROPERTY_DEFAULTS_HH
#define DUMUX_1PNCMIN_PROPERTY_DEFAULTS_HH

#include "model.hh"
#include "indices.hh"
#include "volumevariables.hh"
#include "properties.hh"

#include <dumux/porousmediumflow/compositional/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/porousmediumflow/1pnc/implicit/newtoncontroller.hh>
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/material/spatialparams/implicit1p.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>

namespace Dumux
{

namespace Properties {
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////

/*!
 * \brief Set the property for the number of secondary components.
 * Secondary components are components calculated from
 * primary components by equilibrium relations and
 * do not have mass balance equation on their own.
 * These components are important in the context of bio-mineralization applications.
 * We just forward the number from the fluid system
 *
 */
SET_PROP(OnePNCMin, NumSComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
public:
    static const int value = FluidSystem::numSComponents;
};
/*!
 * \brief Set the property for the number of solid phases, excluding the non-reactive matrix.
 *
 * We just forward the number from the fluid system
 *
 */
SET_PROP(OnePNCMin, NumSPhases)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numSPhases;
};

//physical processes to be considered by the isothermal model
SET_BOOL_PROP(OnePNCMin, EnableAdvection, true);
SET_BOOL_PROP(OnePNCMin, EnableMolecularDiffusion, true);
SET_BOOL_PROP(OnePNCMin, EnableEnergyBalance, false);

/*!
 * \brief Set the property for the number of equations.
 * For each component and each precipitated mineral/solid phase one equation has to
 * be solved.
 */
SET_PROP(OnePNCMin, NumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
//     static const int value = FluidSystem::numComponents + FluidSystem::numSPhases;
    static const int value = FluidSystem::numComponents + FluidSystem::numSComponents; //steamaircao2h2 has 2 components in the fluidphase
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(OnePNCMin, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef CompositionalFluidState<Scalar, FluidSystem> type;
};

//! Use the 2pncmin local residual operator
SET_TYPE_PROP(OnePNCMin,
              LocalResidual,
              OnePNCMinLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(OnePNCMin, Model, OnePNCMinModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(OnePNCMin, VolumeVariables, OnePNCMinVolumeVariables<TypeTag>);

//! the FluxVariables property
// SET_TYPE_PROP(OnePNCMin, FluxVariables, OnePNCMinFluxVariables<TypeTag>);

//! The indices required by the isothermal 2pNcMin model
SET_TYPE_PROP(OnePNCMin, Indices, OnePNCMinIndices <TypeTag, /*PVOffset=*/0>);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(OnePNCMin, SpatialParamsForchCoeff, 0.55);

//set isothermal VolumeVariables
SET_TYPE_PROP(OnePNCMin, IsothermalVolumeVariables, OnePNCMinVolumeVariables<TypeTag>);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(OnePNCMinNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivityAverage<Scalar> type;
};


SET_BOOL_PROP(OnePNCMinNI, NiOutputLevel, 0);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(OnePNCMinNI, IsothermalModel, OnePNCMinModel<TypeTag>);

// set isothermal FluxVariables
//SET_TYPE_PROP(OnePNCMinNI, IsothermalFluxVariables, OnePNCMinFluxVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(OnePNCMinNI, IsothermalVolumeVariables, OnePNCMinVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(OnePNCMinNI, IsothermalLocalResidual, OnePNCMinLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(OnePNCMinNI, IsothermalIndices, OnePNCMinIndices<TypeTag, /*PVOffset=*/0>);

//set isothermal NumEq
SET_PROP(OnePNCMinNI, IsothermalNumEq)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents +FluidSystem::numSComponents;// in NonIsothermal 1 is substracted by default
};
}
}

#endif
