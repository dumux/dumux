// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \file
 * \ingroup SequentialTwoPModel
 * \brief Specifies the properties for immiscible 2p diffusion/pressure models
 */
#ifndef DUMUX_DIFFUSION_PROPERTIES_2P_HH
#define DUMUX_DIFFUSION_PROPERTIES_2P_HH

#include <dumux/porousmediumflow/sequential/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/properties.hh>

namespace Dumux {
namespace Properties {
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(PressureTwoP, INHERITS_FROM(Pressure, SequentialTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
} // end namespace Properties
} // end namespace Dumux

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
SET_TYPE_PROP(PressureTwoP, Model, typename GET_PROP_TYPE(TypeTag, PressureModel));
} // end namespace Properties
} // end namespace Dumux

#endif
