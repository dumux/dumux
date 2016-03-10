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
 * \ingroup FV2p
 * \ingroup Properties
 */
/*!
 * \file
 *
 * \brief Defines the properties required for finite volume pressure models
 */

#ifndef DUMUX_FVPRESSUREVELOCITYPORPERTIES2P_SEQUENTIAL_HH
#define DUMUX_FVPRESSUREVELOCITYPORPERTIES2P_SEQUENTIAL_HH

//Dumux-includes
#include <dumux/porousmediumflow/2p/sequential/diffusion/properties.hh>

namespace Dumux
{

////////////////////////////////
// forward declarations
////////////////////////////////


////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the one-phase problems using a standard finite volume model
NEW_TYPE_TAG(FVPressureVelocityTwoP, INHERITS_FROM(PressureTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

}
}

#include "velocity.hh"
#include "pressurevelocity.hh"

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
//! Set velocity reconstruction implementation standard cell centered finite volume schemes as default
SET_TYPE_PROP( FVPressureVelocityTwoP, Velocity, Dumux::FVVelocity2P<TypeTag> );
//! Set finite volume implementation of the one-phase pressure equation as default pressure model
SET_TYPE_PROP(FVPressureVelocityTwoP, PressureModel, Dumux::FVPressureVelocity2P<TypeTag>);
//! Allow assembling algorithm for the pressure matrix to assemble only from one side of a cell-cell interface
SET_BOOL_PROP(FVPressureVelocityTwoP, VisitFacesOnlyOnce, true);
// \}
}

}

#endif
