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
#ifndef DUMUX_DIFFUSION_PROPERTIES_2P_HH
#define DUMUX_DIFFUSION_PROPERTIES_2P_HH

#include <dumux/decoupled/common/pressureproperties.hh>
#include <dumux/decoupled/2p/2pproperties.hh>

/*!
 * \ingroup Pressure2p
 * \ingroup IMPETProperties
 */
/*!
 * \file
 * \brief Specifies the properties for immiscible 2p diffusion/pressure models
 */
namespace Dumux
{
namespace Properties
{
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(PressureTwoP, INHERITS_FROM(Pressure, DecoupledTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
}
}

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
SET_TYPE_PROP(PressureTwoP, Model, typename GET_PROP_TYPE(TypeTag, PressureModel));
}
}

#endif
