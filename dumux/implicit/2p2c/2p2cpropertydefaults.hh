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

#include "2p2cmodel.hh"
#include "2p2cindices.hh"
#include "2p2cfluxvariables.hh"
#include "2p2cvolumevariables.hh"
#include "2p2cproperties.hh"
#include "2p2clocalresidual.hh"
#include "2p2cnewtoncontroller.hh"

#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/implicit/common/implicitdarcyfluxvariables.hh>
#include <dumux/material/spatialparams/implicitspatialparams.hh>

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

SET_INT_PROP(TwoPTwoC, NumEq, 2); //!< set the number of equations to 2

//! Set the default formulation to pw-sn
SET_INT_PROP(TwoPTwoC,
             Formulation,
             TwoPTwoCFormulation::pwsn);

//! set as default that no component mass balance is replaced by the total mass balance
SET_INT_PROP(TwoPTwoC, ReplaceCompEqIdx, 2);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_PROP(TwoPTwoC, MaterialLawParams)
{
 private:
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;

 public:
    typedef typename MaterialLaw::Params type;
};

//! Use the 2p2c local Jacobian operator for the 2p2c model
SET_TYPE_PROP(TwoPTwoC,
              LocalResidual,
              TwoPTwoCLocalResidual<TypeTag>);

//! Use the 2p2c specific Newton controller for the 2p2c model
SET_TYPE_PROP(TwoPTwoC, NewtonController, TwoPTwoCNewtonController<TypeTag>);

//! the Model property
SET_TYPE_PROP(TwoPTwoC, Model, TwoPTwoCModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(TwoPTwoC, VolumeVariables, TwoPTwoCVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(TwoPTwoC, FluxVariables, TwoPTwoCFluxVariables<TypeTag>);

//! define the base flux variables to realize Darcy flow
SET_TYPE_PROP(TwoPTwoC, BaseFluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//! the upwind weight for the mass conservation equations.
SET_SCALAR_PROP(TwoPTwoC, ImplicitMassUpwindWeight, 1.0);

//! set default mobility upwind weight to 1.0, i.e. fully upwind
SET_SCALAR_PROP(TwoPTwoC, ImplicitMobilityUpwindWeight, 1.0);

//! The indices required by the isothermal 2p2c model
SET_PROP(TwoPTwoC, Indices)
{ private:
    enum { Formulation = GET_PROP_VALUE(TypeTag, Formulation) };
 public:
    typedef TwoPTwoCIndices<TypeTag, Formulation, 0> type;
};

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(TwoPTwoC, SpatialParams, ImplicitSpatialParams<TypeTag>);

//! The model after Millington (1961) is used for the effective diffusivity
SET_PROP(TwoPTwoC, EffectiveDiffusivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef DiffusivityMillingtonQuirk<Scalar> type;
};

// disable velocity output by default
SET_BOOL_PROP(TwoPTwoC, VtkAddVelocity, false);

// enable gravity by default
SET_BOOL_PROP(TwoPTwoC, ProblemEnableGravity, true);

SET_BOOL_PROP(TwoPTwoC, UseMoles, true); //!< Define that mole fractions are used in the balance equations per default


//! default value for the Forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(BoxModel, SpatialParamsForchCoeff, 0.55);}

}

#endif
