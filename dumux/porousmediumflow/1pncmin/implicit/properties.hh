// -**- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \ingroup OnePNCMinModel
 *
 * \file
 *
 * \brief Defines the properties required for the one-phase n-component mineralization
 *        fully implicit model.
 */
#ifndef DUMUX_1PNCMIN_PROPERTIES_HH
#define DUMUX_1PNCMIN_PROPERTIES_HH

#include <dumux/porousmediumflow/1pnc/model.hh>

namespace Dumux {
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the isothermal two phase n component mineralisation problems
NEW_TYPE_TAG(OnePNCMin, INHERITS_FROM(OnePNC));
NEW_TYPE_TAG(BoxOnePNCMin, INHERITS_FROM(BoxModel, OnePNCMin));
NEW_TYPE_TAG(CCOnePNCMin, INHERITS_FROM(CCModel, OnePNCMin));

//! The type tags for the corresponding non-isothermal problems
NEW_TYPE_TAG(OnePNCMinNI, INHERITS_FROM(OnePNCMin, NonIsothermal));
NEW_TYPE_TAG(BoxOnePNCMinNI, INHERITS_FROM(BoxModel, OnePNCMinNI));
NEW_TYPE_TAG(CCOnePNCMinNI, INHERITS_FROM(CCModel, OnePNCMinNI));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumSPhases); //!< Number of solid phases in the system
NEW_PROP_TAG(NumFSPhases); //!< Number of fluid and solid phases in the system
NEW_PROP_TAG(NumSComponents); //!< Number of solid components in the system
NEW_PROP_TAG(NumPSComponents); //!< Number of fluid and solid components in the system
NEW_PROP_TAG(NumTraceComponents); //!< Number of trace fluid components which are not considered in the calculation of the phase density
NEW_PROP_TAG(NumSecComponents); //!< Number of secondary components which are not primary variables
NEW_PROP_TAG(OnePNCMinIndices); //!< Enumerations for the 2pncMin models
NEW_PROP_TAG(SpatialParamsForchCoeff); //!< Property for the forchheimer coefficient

}
}

#endif
