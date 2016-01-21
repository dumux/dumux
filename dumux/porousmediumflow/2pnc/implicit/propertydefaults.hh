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
#include "fluxvariables.hh"
#include "volumevariables.hh"
#include "properties.hh"
#include "newtoncontroller.hh"

#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
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
//! The major components belonging to the existing phases are mentioned here e.g., 2 for water and air being the major component in the liquid and gas phases in a 2 phase system
SET_PROP(TwoPNC, NumMajorComponents)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "The model is restricted to two-phases, thus number of major components must also be two.");
};

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
SET_PROP(TwoPNC, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

//! Set the default formulation to pl-Sg: This can be over written in the problem.
SET_INT_PROP(TwoPNC, Formulation, TwoPNCFormulation::plSg);

//! Set the property for the material parameters by extracting it from the material law.
SET_PROP(TwoPNC, MaterialLawParams)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(MaterialLaw)) MaterialLaw;

public:
    typedef typename MaterialLaw::Params type;
};

//! Use the 2pnc local residual
SET_TYPE_PROP(TwoPNC, LocalResidual, TwoPNCLocalResidual<TypeTag>);

//! Use the 2pnc newton controller
SET_TYPE_PROP(TwoPNC, NewtonController, TwoPNCNewtonController<TypeTag>);

//! the Model property
SET_TYPE_PROP(TwoPNC, Model, TwoPNCModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(TwoPNC, VolumeVariables, TwoPNCVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(TwoPNC, FluxVariables, TwoPNCFluxVariables<TypeTag>);

//! define the base flux variables to realize Darcy flow
SET_TYPE_PROP(TwoPNC, BaseFluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//! the upwind weight for the mass conservation equations.
SET_SCALAR_PROP(TwoPNC, ImplicitMassUpwindWeight, 1.0);

//! Set default mobility upwind weight to 1.0, i.e. fully upwind
SET_SCALAR_PROP(TwoPNC, ImplicitMobilityUpwindWeight, 1.0);

//! The indices required by the isothermal 2pnc model
SET_TYPE_PROP(TwoPNC, Indices, TwoPNCIndices <TypeTag, /*PVOffset=*/0>);

//! Use the ImplicitSpatialParams by default
SET_TYPE_PROP(TwoPNC, SpatialParams, ImplicitSpatialParams<TypeTag>);

//! Enable gravity by default
SET_BOOL_PROP(TwoPNC, ProblemEnableGravity, true);

//! Disable velocity output by default
SET_BOOL_PROP(TwoPNC, VtkAddVelocity, false);

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

// set isothermal FluxVariables
SET_TYPE_PROP(TwoPNCNI, IsothermalFluxVariables, TwoPNCFluxVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(TwoPNCNI, IsothermalVolumeVariables, TwoPNCVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(TwoPNCNI, IsothermalLocalResidual, TwoPNCLocalResidual<TypeTag>);

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
