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
 * \ingroup BoxProperties
 * \ingroup OnePBoxModel
 * \file
 *
 * \brief Defines the properties required for the one-phase BOX model.
 */
#ifndef DUMUX_1P_PROPERTY_DEFAULTS_HH
#define DUMUX_1P_PROPERTY_DEFAULTS_HH

#include <dumux/implicit/box/boxproperties.hh>

#include "1pmodel.hh"
#include "1plocalresidual.hh"
#include "1pvolumevariables.hh"
#include <dumux/implicit/box/boxdarcyfluxvariables.hh>
#include "1pindices.hh"

#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/1pfluidsystem.hh>
#include <dumux/material/spatialparams/boxspatialparams1p.hh>

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
SET_TYPE_PROP(OneP,
              LocalResidual,
              OnePLocalResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(OneP, Model, OnePBoxModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(OneP, VolumeVariables, OnePVolumeVariables<TypeTag>);

//! the FluxVariables property
SET_TYPE_PROP(OneP, FluxVariables, BoxDarcyFluxVariables<TypeTag>);

//! The indices required by the isothermal single-phase model
SET_TYPE_PROP(OneP, Indices, OnePIndices);

//! The spatial parameters to be employed. 
//! Use BoxSpatialParamsOneP by default.
SET_TYPE_PROP(OneP, SpatialParams, BoxSpatialParamsOneP<TypeTag>);

//! The weight of the upwind control volume when calculating
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(OneP, ImplicitMassUpwindWeight, 0.5);

//! weight for the upwind mobility in the velocity calculation
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(OneP, ImplicitMobilityUpwindWeight, 0.5);

//! The fluid system to use by default
SET_TYPE_PROP(OneP, FluidSystem, Dumux::FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

SET_PROP(OneP, Fluid)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

// enable gravity by default
SET_BOOL_PROP(OneP, ProblemEnableGravity, true);

//! default value for the forchheimer coefficient
// Source: Ward, J.C. 1964 Turbulent flow in porous media. ASCE J. Hydraul. Div 90.
//        Actually the Forchheimer coefficient is also a function of the dimensions of the
//        porous medium. Taking it as a constant is only a first approximation
//        (Nield, Bejan, Convection in porous media, 2006, p. 10)
SET_SCALAR_PROP(BoxModel, SpatialParamsForchCoeff, 0.55);

// \}
} // end namepace Properties

} // end namepace Dumux

#endif
