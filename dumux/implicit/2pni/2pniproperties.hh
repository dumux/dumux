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
 *  \ingroup TwoPNIModel
 * \file
 *
 * \brief Defines the properties required for the non-isotherm two-phase fully implicit model.
 */
#ifndef DUMUX_2PNI_PROPERTIES_HH
#define DUMUX_2PNI_PROPERTIES_HH

#warning You are including the old 2pni model which will be removed after 2.6. See CHANGELOG for details.

#include <dumux/implicit/2p/2pproperties.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit non-isothermal two-phase problems
NEW_TYPE_TAG(TwoPNI, INHERITS_FROM(TwoP));
NEW_TYPE_TAG(BoxTwoPNI, INHERITS_FROM(BoxModel, TwoPNI));
NEW_TYPE_TAG(CCTwoPNI, INHERITS_FROM(CCModel, TwoPNI));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(ThermalConductivityModel);   //!< The model for the effective thermal conductivity

}
}

#endif
