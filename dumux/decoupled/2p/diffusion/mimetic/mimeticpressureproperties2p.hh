// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \ingroup MimeticPressure2p
 * \ingroup Properties
 */
/*!
 * \file
 *
 * \brief Defines the properties required for a two-phase mimetic finite differences model.
 */

#ifndef DUMUX_MIMETICPROPERTIES2P_DECOUPLED_HH
#define DUMUX_MIMETICPROPERTIES2P_DECOUPLED_HH

//Dumux-includes
#include <dumux/decoupled/2p/diffusion/diffusionproperties2p.hh>
#include <dumux/decoupled/common/mimetic/mimeticproperties.hh>
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

//! The type tag for two-phase problems using a mimetic finite differences method.
NEW_TYPE_TAG(MimeticPressureTwoP, INHERITS_FROM(PressureTwoP, Mimetic))
;

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
}
}

#include <dumux/decoupled/2p/diffusion/mimetic/mimeticpressure2p.hh>
#include <dumux/decoupled/2p/diffusion/mimetic/mimeticgroundwater.hh>

namespace Dumux
{
namespace Properties
{
//! Set mimetic finite differences implementation of the two-phase pressure equation as default pressure model
SET_TYPE_PROP(MimeticPressureTwoP, PressureModel, MimeticPressure2P<TypeTag>);
//! Set the local stiffness implementation for the two-phase model
SET_TYPE_PROP(MimeticPressureTwoP, LocalStiffness, MimeticGroundwaterEquationLocalStiffness<TypeTag>);
}
}

#endif
