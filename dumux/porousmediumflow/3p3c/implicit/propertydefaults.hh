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
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        three-phase three-component fully implicit model.
 */
#ifndef DUMUX_3P3C_PROPERTY_DEFAULTS_HH
#define DUMUX_3P3C_PROPERTY_DEFAULTS_HH

#include "indices.hh"

#include "model.hh"
#include "volumevariables.hh"
#include "properties.hh"
#include "newtoncontroller.hh"
#include "localresidual.hh"
#include "primaryvariableswitch.hh"

#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/3p/thermalconductivitysomerton3p.hh>

namespace Dumux
{

namespace Properties {
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

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(ThreePThreeC, MaterialLawParams, typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

//! The local residual function of the conservation equations
SET_TYPE_PROP(ThreePThreeC, LocalResidual, ThreePThreeCLocalResidual<TypeTag>);

//! Enable advection
SET_BOOL_PROP(ThreePThreeC, EnableAdvection, true);

//! Enable molecular diffusion
SET_BOOL_PROP(ThreePThreeC, EnableMolecularDiffusion, true);

//! Isothermal model by default
SET_BOOL_PROP(ThreePThreeC, EnableEnergyBalance, false);

//! Use the 3p3c specific newton controller for the 3p3c model
SET_TYPE_PROP(ThreePThreeC, NewtonController, ThreePThreeCNewtonController<TypeTag>);

//! the Model property
SET_TYPE_PROP(ThreePThreeC, Model, ThreePThreeCModel<TypeTag>);

//! The primary variable switch for the 3p3c model
SET_TYPE_PROP(ThreePThreeC, PrimaryVariableSwitch, ThreePThreeCPrimaryVariableSwitch<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(ThreePThreeC, VolumeVariables, ThreePThreeCVolumeVariables<TypeTag>);

//! the upwind factor for the mobility.
SET_SCALAR_PROP(ThreePThreeC, ImplicitMassUpwindWeight, 1.0);

//! set default mobility upwind weight to 1.0, i.e. fully upwind
SET_SCALAR_PROP(ThreePThreeC, ImplicitMobilityUpwindWeight, 1.0);

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

// disable velocity output by default
SET_BOOL_PROP(ThreePThreeC, VtkAddVelocity, false);

// enable gravity by default
SET_BOOL_PROP(ThreePThreeC, ProblemEnableGravity, true);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(BoxModel, SpatialParamsForchCoeff, 0.55);

/*!
 * \brief default value for tortuosity value (tau) used in macroscopic diffusion
 *
 * Value is 0.5 according to Carman 1937: <i>Fluid flow through granular beds</i>
 * \cite carman1937
 */
SET_SCALAR_PROP(ThreePThreeC, TauTortuosity, 0.5);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(ThreePThreeCNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

//! temperature is already written by the isothermal model
SET_BOOL_PROP(ThreePThreeCNI, NiOutputLevel, 0);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(ThreePThreeCNI, IsothermalModel, ThreePThreeCModel<TypeTag>);

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

}

}

#endif
