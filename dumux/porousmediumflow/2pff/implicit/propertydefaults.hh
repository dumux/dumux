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
#ifndef DUMUX_2P_FRACFLOW_PROPERTY_DEFAULTS_HH
#define DUMUX_2P_FRACFLOW_PROPERTY_DEFAULTS_HH

#include "properties.hh"
#include "localresidual.hh"
#include "volumevariables.hh"
#include "darcyslaw.hh"
#include "upwindscheme.hh"
#include <dumux/porousmediumflow/2p/implicit/indices.hh>

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property defaults
//////////////////////////////////////////////////////////////////
//! Set the default formulation to pWsN
SET_INT_PROP(TwoPFractionalFlow,
             Formulation,
             TwoPFormulation::pnsw);

//! We have a special local resdiual
SET_TYPE_PROP(TwoPFractionalFlow,
              LocalResidual,
              TwoPFractionalFlowLocalResidual<TypeTag>);

//! We have a special advection type
SET_TYPE_PROP(TwoPFractionalFlow,
              AdvectionType,
              FractionalFlowDarcysLaw<TypeTag>);

//! We have special volume variables (computing pressures / saturations from the total velocity depending on the privar?)
SET_TYPE_PROP(TwoPFractionalFlow,
              VolumeVariables,
              TwoPFractionalFlowVolumeVariables<TypeTag>);

//! We have a special upwind scheme
SET_TYPE_PROP(TwoPFractionalFlow,
              UpwindScheme,
              TwoPFractionalFlowUpwindScheme<TypeTag>);

} // end namespace Properties

} // end namespace Dumux

#endif
