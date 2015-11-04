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
 * \ingroup TwoPModel
 * \ingroup ImplicitProperties
 * \ingroup Properties
 * \file
 *
 * \brief Defines default values for the properties required by the
 *        two-phase fully implicit model.
 */
#ifndef DUMUX_2P_PROPERTY_DEFAULTS_HH
#define DUMUX_2P_PROPERTY_DEFAULTS_HH

#include "2pproperties.hh"
#include "2pmodel.hh"
#include "2pindices.hh"
#include "2pvolumevariables.hh"
#include "2plocalresidual.hh"

#include <dumux/implicit/nonisothermal/nipropertydefaults.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/2pimmisciblefluidsystem.hh>
#include <dumux/material/fluidstates/immisciblefluidstate.hh>
#include <dumux/material/spatialparams/implicitspatialparams.hh>
#include <dumux/implicit/common/implicitdarcyfluxvariables.hh>
#include <dumux/material/fluidmatrixinteractions/2p/thermalconductivitysomerton.hh>

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
//! Set the number of equations to 2
SET_INT_PROP(TwoP, NumEq, 2);

//! The number of phases in the 2p model is 2
SET_INT_PROP(TwoP, NumPhases, 2);

//! Set the default formulation to pWsN
SET_INT_PROP(TwoP,
             Formulation,
             TwoPFormulation::pwsn);

//! Use the 2p local jacobian operator for the 2p model
SET_TYPE_PROP(TwoP,
              LocalResidual,
              TwoPLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(TwoP, Model, TwoPModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(TwoP, VolumeVariables, TwoPVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(TwoP, FluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//! the upwind weight for the mass conservation equations.
SET_SCALAR_PROP(TwoP, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(TwoP, ImplicitMobilityUpwindWeight, 1.0);

//! The indices required by the isothermal 2p model
SET_TYPE_PROP(TwoP,
              Indices,
              TwoPIndices<TypeTag, GET_PROP_VALUE(TypeTag, Formulation), 0>);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParams by default.
SET_TYPE_PROP(TwoP, SpatialParams, ImplicitSpatialParams<TypeTag>);

/*!
 * \brief Set the property for the material parameters by extracting
 *        it from the material law.
 */
SET_TYPE_PROP(TwoP,
              MaterialLawParams,
              typename GET_PROP_TYPE(TypeTag, MaterialLaw)::Params);

SET_PROP(TwoP, WettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(TwoP, NonwettingPhase)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

SET_PROP(TwoP, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, WettingPhase) WettingPhase;
    typedef typename GET_PROP_TYPE(TypeTag, NonwettingPhase) NonwettingPhase;

public:
    typedef Dumux::FluidSystems::TwoPImmiscible<Scalar,
                                                WettingPhase,
                                                NonwettingPhase> type;
};

SET_PROP(TwoP, FluidState)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
public:
    typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};

// disable velocity output by default
SET_BOOL_PROP(TwoP, VtkAddVelocity, false);

// enable gravity by default
SET_BOOL_PROP(TwoP, ProblemEnableGravity, true);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(TwoP, SpatialParamsForchCoeff, 0.55);

//! Somerton is used as default model to compute the effective thermal heat conductivity
SET_PROP(TwoPNI, ThermalConductivityModel)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
public:
    typedef ThermalConductivitySomerton<Scalar, Indices> type;
};

//! temperature is already written by the isothermal model
SET_BOOL_PROP(TwoPNI, NiOutputLevel, 0);

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(TwoPNI, IsothermalModel, TwoPModel<TypeTag>);

// set isothermal FluxVariables
SET_TYPE_PROP(TwoPNI, IsothermalFluxVariables, ImplicitDarcyFluxVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(TwoPNI, IsothermalVolumeVariables, TwoPVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(TwoPNI, IsothermalLocalResidual, TwoPLocalResidual<TypeTag>);

//set isothermal Indices
SET_PROP(TwoPNI, IsothermalIndices)
{
private:
    enum { Formulation = GET_PROP_VALUE(TypeTag, Formulation) };
public:
    typedef TwoPIndices<TypeTag, Formulation, 0> type;
};

//set isothermal NumEq
SET_INT_PROP(TwoPNI, IsothermalNumEq, 2);

}


}

#endif
