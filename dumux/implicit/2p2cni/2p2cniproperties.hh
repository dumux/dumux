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
 * \ingroup TwoPTwoCNIModel
 * \file
 *
 * \brief Defines the properties required for the non-isothermal two-phase,
 * two-component fully implicit model.
 */
#ifndef DUMUX_2P2CNI_PROPERTIES_HH
#define DUMUX_2P2CNI_PROPERTIES_HH

#include <dumux/implicit/2p2c/2p2cproperties.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit non-isothermal two-phase two-component problems
NEW_TYPE_TAG(TwoPTwoCNI, INHERITS_FROM(TwoPTwoC));
NEW_TYPE_TAG(BoxTwoPTwoCNI, INHERITS_FROM(BoxModel, TwoPTwoCNI));
NEW_TYPE_TAG(CCTwoPTwoCNI, INHERITS_FROM(CCModel, TwoPTwoCNI));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(ThermalConductivityModel);   //!< The model for the effective thermal conductivity
}
}

#endif
