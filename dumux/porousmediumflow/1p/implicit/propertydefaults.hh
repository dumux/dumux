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
 * \ingroup OnePModel
 * \file
 *
 * \brief Defines the properties required for the one-phase fully implicit model.
 */
#ifndef DUMUX_1P_PROPERTY_DEFAULTS_HH
#define DUMUX_1P_PROPERTY_DEFAULTS_HH

#include "properties.hh"

#include "model.hh"
#include "volumevariables.hh"
#include "indices.hh"

#include <dumux/porousmediumflow/immiscible/localresidual.hh>
#include <dumux/porousmediumflow/nonisothermal/implicit/propertydefaults.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/1p.hh>
#include <dumux/material/spatialparams/implicit1p.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>

namespace Dumux
{
// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {
SET_INT_PROP(OneP, NumEq, 1); //!< set the number of equations to 1
SET_INT_PROP(OneP, NumPhases, 1); //!< The number of phases in the 1p model is 1

//! The local residual function
SET_TYPE_PROP(OneP, LocalResidual, ImmiscibleLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(OneP, Model, OnePModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(OneP, VolumeVariables, OnePVolumeVariables<TypeTag>);

//! Enable advection
SET_BOOL_PROP(OneP, EnableAdvection, true);

//! The one-phase model has no molecular diffusion
SET_BOOL_PROP(OneP, EnableMolecularDiffusion, false);

//! Isothermal model by default
SET_BOOL_PROP(OneP, EnableEnergyBalance, false);

//! The indices required by the isothermal single-phase model
SET_TYPE_PROP(OneP, Indices, OnePIndices);

//! The spatial parameters to be employed.
//! Use ImplicitSpatialParamsOneP by default.
SET_TYPE_PROP(OneP, SpatialParams, ImplicitSpatialParamsOneP<TypeTag>);

//! The weight of the upwind control volume when calculating
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(OneP, ImplicitMassUpwindWeight, 1.0);

//! weight for the upwind mobility in the velocity calculation
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(OneP, ImplicitMobilityUpwindWeight, 1.0);

//! The fluid system to use by default
SET_TYPE_PROP(OneP, FluidSystem, FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

SET_PROP(OneP, Fluid)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, NullComponent<Scalar> > type;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(OneP, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef ImmiscibleFluidState<Scalar, FluidSystem> type;
};

// disable velocity output by default

// enable gravity by default
SET_BOOL_PROP(OneP, ProblemEnableGravity, true);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(OneP, SpatialParamsForchCoeff, 0.55);

//! average is used as default model to compute the effective thermal heat conductivity
SET_PROP(OnePNI, ThermalConductivityModel)
{ private :
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
  public:
    typedef ThermalConductivityAverage<Scalar> type;
};

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(OnePNI, IsothermalModel, OnePModel<TypeTag>);

//set isothermal VolumeVariables
SET_TYPE_PROP(OnePNI, IsothermalVolumeVariables, OnePVolumeVariables<TypeTag>);

//set isothermal LocalResidual
SET_TYPE_PROP(OnePNI, IsothermalLocalResidual, ImmiscibleLocalResidual<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(OnePNI, IsothermalIndices, OnePIndices);

//set isothermal NumEq
SET_INT_PROP(OnePNI, IsothermalNumEq, 1);

// \}
} // end namespace Properties

} // end namespace Dumux

#endif
