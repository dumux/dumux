// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2011 by Markus Wolff                                      *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

/*!
 * \ingroup IMPES
 * \ingroup Properties
 */
/*!
 * \file
 *
 * \brief Defines the properties required for finite volume pressure models in a twophase sequential model.
 */

#ifndef DUMUX_FVPRESSUREPORPERTIES2P_ADAPTIVE_DECOUPLED_HH
#define DUMUX_FVPRESSUREPORPERTIES2P_ADAPTIVE_DECOUPLED_HH

//Dumux-includes
#include <dumux/decoupled/2p/diffusion/diffusionproperties2p.hh>

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

//! The type tag for the two-phase problems
NEW_TYPE_TAG(FVPressureTwoPAdaptive, INHERITS_FROM(PressureTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

}
}

#include "fvvelocity2padaptive.hh"
#include "fvpressure2padaptive.hh"

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
SET_TYPE_PROP( FVPressureTwoPAdaptive, Velocity, Dumux::FVVelocity2PAdaptive<TypeTag> );
SET_TYPE_PROP(FVPressureTwoPAdaptive, PressureModel, Dumux::FVPressure2PAdaptive<TypeTag>);
SET_BOOL_PROP(FVPressureTwoPAdaptive, VisitFacesOnlyOnce, true);
// \}
}

}

#endif
