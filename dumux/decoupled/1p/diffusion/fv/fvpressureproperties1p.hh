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

#ifndef DUMUX_FVPRESSUREPORPERTIES1P_DECOUPLED_HH
#define DUMUX_FVPRESSUREPORPERTIES1P_DECOUPLED_HH

//Dumux-includes
#include <dumux/decoupled/1p/diffusion/diffusionproperties1p.hh>

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
NEW_TYPE_TAG(FVPressureOneP, INHERITS_FROM(PressureOneP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

}
}

#include "fvvelocity1p.hh"
#include "fvpressure1p.hh"

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
SET_TYPE_PROP( FVPressureOneP, Velocity, Dumux::FVVelocity1P<TypeTag> );
SET_TYPE_PROP(FVPressureOneP, PressureModel, Dumux::FVPressure1P<TypeTag>);
//! Faces are only regarded from one side and not from both cells
SET_BOOL_PROP(FVPressureOneP, VisitFacesOnlyOnce, true);
// \}
}

}

#endif
