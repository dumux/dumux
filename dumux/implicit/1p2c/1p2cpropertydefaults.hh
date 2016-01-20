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
 * \ingroup OnePTwoCModel
 * \file
 *
 * \brief Defines some default values for the properties of the the
 *        single-phase, two-component fully implicit model.
 */

#ifndef DUMUX_1P2C_PROPERTY_DEFAULTS_HH
#define DUMUX_1P2C_PROPERTY_DEFAULTS_HH

#include "1p2cproperties.hh"
#include "1p2cmodel.hh"
#include "1p2clocalresidual.hh"
#include "1p2cvolumevariables.hh"
#include "1p2cfluxvariables.hh"
#include "1p2cindices.hh"

#include <dumux/porousmediumflow/implicit/nonisothermal/propertydefaults.hh>
#include <dumux/material/spatialparams/implicitspatialparams1p.hh>
#include <dumux/material/fluidmatrixinteractions/diffusivitymillingtonquirk.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/fluidstates/compositionalfluidstate.hh>

namespace Dumux
{
// \{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property values
//////////////////////////////////////////////////////////////////


SET_INT_PROP(OnePTwoC, NumEq, 2); //!< set the number of equations to 2
SET_INT_PROP(OnePTwoC, NumPhases, 1); //!< The number of phases in the 1p2c model is 1
SET_INT_PROP(OnePTwoC, NumComponents, 2); //!< The number of components in the 1p2c model is 2
SET_SCALAR_PROP(OnePTwoC, Scaling, 1); //!< Scaling of the model is set to 1 by default
SET_BOOL_PROP(OnePTwoC, UseMoles, true); //!< Define that mole fractions are used in the balance equations

//! Use the 1p2c local residual function for the 1p2c model
SET_TYPE_PROP(OnePTwoC, LocalResidual, OnePTwoCLocalResidual<TypeTag>);

//! define the model
SET_TYPE_PROP(OnePTwoC, Model, OnePTwoCModel<TypeTag>);

//! define the VolumeVariables
SET_TYPE_PROP(OnePTwoC, VolumeVariables, OnePTwoCVolumeVariables<TypeTag>);

//! define the FluxVariables
SET_TYPE_PROP(OnePTwoC, FluxVariables, OnePTwoCFluxVariables<TypeTag>);

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(OnePTwoC, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::CompositionalFluidState<Scalar, FluidSystem> type;
};

//! set default upwind weight to 1.0, i.e. fully upwind
SET_SCALAR_PROP(OnePTwoC, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
SET_SCALAR_PROP(OnePTwoC, ImplicitMobilityUpwindWeight, 1.0);

//! Set the indices used by the 1p2c model
SET_TYPE_PROP(OnePTwoC, Indices, OnePTwoCIndices<TypeTag>);
//! The spatial parameters to be employed.
//! Use ImplicitSpatialParamsOneP by default.
SET_TYPE_PROP(OnePTwoC, SpatialParams, ImplicitSpatialParamsOneP<TypeTag>);

//! The model after Millington (1961) is used for the effective diffusivity
SET_PROP(OnePTwoC, EffectiveDiffusivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
 public:
    typedef DiffusivityMillingtonQuirk<Scalar> type;
};

//! Set the phaseIndex per default to zero (important for two-phase fluidsystems).
SET_INT_PROP(OnePTwoC, PhaseIdx, 0);

// disable velocity output by default
SET_BOOL_PROP(OnePTwoC, VtkAddVelocity, false);

// enable gravity by default
SET_BOOL_PROP(OnePTwoC, ProblemEnableGravity, true);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(OnePTwoC, SpatialParamsForchCoeff, 0.55);

//! average is used as default model to compute the effective thermal heat conductivity
SET_PROP(OnePTwoCNI, ThermalConductivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
  public:
    typedef ThermalConductivityAverage<Scalar> type;
};

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(OnePTwoCNI, IsothermalModel, OnePTwoCModel<TypeTag>);

// set isothermal FluxVariables
SET_TYPE_PROP(OnePTwoCNI, IsothermalFluxVariables, OnePTwoCFluxVariables<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(OnePTwoCNI, IsothermalVolumeVariables, OnePTwoCVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(OnePTwoCNI, IsothermalLocalResidual, OnePTwoCLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(OnePTwoCNI, IsothermalIndices, OnePTwoCIndices<TypeTag>);

//set isothermal NumEq
SET_INT_PROP(OnePTwoCNI, IsothermalNumEq, 2);


}
// \}
}

#endif

