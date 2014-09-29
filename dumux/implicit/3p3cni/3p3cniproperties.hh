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
 * \ingroup ThreePThreeCNIModel
 * \file
 *
 * \brief Defines the properties required for the non-isothermal three-phase
 *        three-component fully implicit model.
 */
#ifndef DUMUX_3P3CNI_PROPERTIES_HH
#define DUMUX_3P3CNI_PROPERTIES_HH

#warning You are including the old 2pni model which will be removed after 2.6. See CHANGELOG for details.

#include <dumux/implicit/3p3c/3p3cproperties.hh>

namespace Dumux
{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the implicit non-isothermal three-phase three-component problems
NEW_TYPE_TAG(ThreePThreeCNI, INHERITS_FROM(ThreePThreeC));
NEW_TYPE_TAG(BoxThreePThreeCNI, INHERITS_FROM(BoxModel, ThreePThreeCNI));
NEW_TYPE_TAG(CCThreePThreeCNI, INHERITS_FROM(CCModel, ThreePThreeCNI));
}
}

#endif
