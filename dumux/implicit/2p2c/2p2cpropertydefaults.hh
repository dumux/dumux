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
 * \file
 *
 * \brief Defines default values for most properties required by the
 *        two-phase two-component fully implicit model.
 */
#ifndef DUMUX_2P2C_PROPERTY_DEFAULTS_HH
#define DUMUX_2P2C_PROPERTY_DEFAULTS_HH

#include "2p2cproperties.hh"
#include "2p2cmodel.hh"
#include "2p2cindices.hh"
#include "2p2cfluxvariables.hh"
#include "2p2cvolumevariables.hh"
#include "2p2clocalresidual.hh"
#include "2p2cnewtoncontroller.hh"

#include <dumux/implicit/nonisothermal/nipropertydefaults.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/porousmediumflow/implicit/darcyfluxvariables.hh>
#include <dumux/material/spatialparams/implicitspatialparams.hh>
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
 * We just forward the number from the fluid system and use a static
 * assert to make sure it is 2.
 */
SET_PROP(TwoPTwoC, NumComponents)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

 public:
    static const int value = FluidSystem::numComponents;

    static_assert(value == 2,
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
 private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

 public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "Only fluid systems with 2 phases are supported by the 2p-2c model!");
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(TwoPTwoC, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

//! Set the number of equations to 2
SET_INT_PROP(TwoPTwoC, NumEq, 2);

//! Set the default formulation to pw-sn
SET_INT_PROP(TwoPTwoC,
             Formulation,
             TwoPTwoCFormulation::pwsn);

//! Set as default that no component mass balance is replaced by the total mass balance
SET_INT_PROP(TwoPTwoC, ReplaceCompEqIdx, 2);

//! Set the property for the material parameters by extracting it from the material law.
SET_PROP(TwoPTwoC, MaterialLawParams)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

 public:
    typedef typename MaterialLaw::Params type;
};

//! Use the 2p2c local residual operator
SET_TYPE_PROP(TwoPTwoC,
              LocalResidual,
              TwoPTwoCLocalResidual<TypeTag>);

//! Use the 2p2c Newton controller
SET_TYPE_PROP(TwoPTwoC, NewtonController, TwoPTwoCNewtonController<TypeTag>);

//! Use the 2p2c model
SET_TYPE_PROP(TwoPTwoC, Model, TwoPTwoCModel<TypeTag>);

//! Use the 2p2c VolumeVariables
SET_TYPE_PROP(TwoPTwoC, VolumeVariables, TwoPTwoCVolumeVariables<TypeTag>);

//! Use the 2p2c FluxVariables
SET_TYPE_PROP(TwoPTwoC, FluxVariables, TwoPTwoCFluxVariables<TypeTag>);

//! Set the BaseFluxVariables to realize Darcy flow
SET_TYPE_PROP(TwoPTwoC, BaseFluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//! Set the upwind weight for the mass conservation equations
SET_SCALAR_PROP(TwoPTwoC, ImplicitMassUpwindWeight, 1.0);

//! Set default mobility upwind weight to 1.0, i.e. fully upwind
SET_SCALAR_PROP(TwoPTwoC, ImplicitMobilityUpwindWeight, 1.0);

//! Set the indices required by the isothermal 2p2c
SET_PROP(TwoPTwoC, Indices)
{ private:
    enum { Formulation = GET_PROP_VALUE(TypeTag, Formulation) };
 public:
    typedef TwoPTwoCIndices<TypeTag, Formulation, 0> type;
};

//! Use the ImplicitSpatialParams by default
SET_TYPE_PROP(TwoPTwoC, SpatialParams, ImplicitSpatialParams<TypeTag>);

//! Use the model after Millington (1961) for the effective diffusivity
SET_PROP(TwoPTwoC, EffectiveDiffusivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef DiffusivityMillingtonQuirk<Scalar> type;
};

//! Disable velocity output by default
SET_BOOL_PROP(TwoPTwoC, VtkAddVelocity, false);

//! Enable gravity by default
SET_BOOL_PROP(TwoPTwoC, ProblemEnableGravity, true);

//! Use mole fractions in the balance equations by default
SET_BOOL_PROP(TwoPTwoC, UseMoles, true);

//! Determines whether the constraint solver is used
SET_BOOL_PROP(TwoPTwoC, UseConstraintSolver, true);

//! Set default value for the Forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(TwoPTwoC, SpatialParamsForchCoeff, 0.55);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(TwoPTwoCNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

//! temperature is already written by the isothermal model
SET_BOOL_PROP(TwoPTwoCNI, NiOutputLevel, 0);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(TwoPTwoCNI, IsothermalModel, TwoPTwoCModel<TypeTag>);

// set isothermal FluxVariables
SET_TYPE_PROP(TwoPTwoCNI, IsothermalFluxVariables, TwoPTwoCFluxVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(TwoPTwoCNI, IsothermalVolumeVariables, TwoPTwoCVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(TwoPTwoCNI, IsothermalLocalResidual, TwoPTwoCLocalResidual<TypeTag>);

//set isothermal Indices
SET_PROP(TwoPTwoCNI, IsothermalIndices)
{
private:
    enum { Formulation = GET_PROP_VALUE(TypeTag, Formulation) };
public:
    typedef TwoPTwoCIndices<TypeTag, Formulation, 0> type;
};

//set isothermal NumEq
SET_INT_PROP(TwoPTwoCNI, IsothermalNumEq, 2);

}

}

#endif
