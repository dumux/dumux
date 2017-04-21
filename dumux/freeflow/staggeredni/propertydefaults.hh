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
#ifndef DUMUX_NAVIER_STOKES_NI_PROPERTY_DEFAULTS_HH
#define DUMUX_NAVIER_STOKES_NI_PROPERTY_DEFAULTS_HH

// TODO clean-up
#include "properties.hh"

#include "model.hh"
#include "volumevariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "fluxvariables.hh"
#include "../staggered/problem.hh"
#include "../staggered/propertydefaults.hh"

#include <dumux/implicit/staggered/localresidual.hh>
//#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
//#include <dumux/material/components/nullcomponent.hh>
//#include <dumux/material/fluidsystems/1p.hh>

#include <dumux/material/fluidstates/immiscible.hh>


namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(FluxVariables);
NEW_PROP_TAG(FluxVariablesCache);
}
// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the non-isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

SET_INT_PROP(NavierStokesNI, NumEqCellCenter, 2);

//! the VolumeVariables property
SET_TYPE_PROP(NavierStokesNI, VolumeVariables, NavierStokesNIVolumeVariables<TypeTag>);
SET_TYPE_PROP(NavierStokesNI, Model, NavierStokesNIModel<TypeTag>);
SET_TYPE_PROP(NavierStokesNI, Indices, NavierStokesNIIndices<TypeTag>);

SET_BOOL_PROP(NavierStokesNI, EnableEnergyBalanceStokes, true);

SET_BOOL_PROP(NavierStokesNI, UseMoles, true);

SET_TYPE_PROP(NavierStokesNI, HeatConductionType, FouriersLaw<TypeTag>);

SET_INT_PROP(NavierStokesNI, PhaseIdx, 0); //!< Defines the phaseIdx

} // end namespace Properties

} // end namespace Dumux

#endif
