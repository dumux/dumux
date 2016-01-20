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
 * \file
 * \brief This file declares and defines the properties required by
 *        the kinetic modules the M-phase N-component model.
 */
#ifndef DUMUX_MPNC_PROPERTIES_KINETIC_HH
#define DUMUX_MPNC_PROPERTIES_KINETIC_HH

#include "properties.hh"

namespace Dumux
{
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
NEW_TYPE_TAG(BoxMPNCKinetic, INHERITS_FROM(BoxMPNC));
//NEW_TYPE_TAG(CCMPNCKinetic, INHERITS_FROM(CCMPNC));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG(AwnSurface);   //!< The material law which ought to be used (extracted from the spatialParams)
NEW_PROP_TAG(AwnSurfaceParams); //!< The context material law (extracted from the spatialParams)
NEW_PROP_TAG(AwsSurface);   //!< The material law which ought to be used (extracted from the spatialParams)
NEW_PROP_TAG(AwsSurfaceParams); //!< The context material law (extracted from the spatialParams)
NEW_PROP_TAG(AnsSurface);   //!< The material law which ought to be used (extracted from the spatialParams)
NEW_PROP_TAG(AnsSurfaceParams); //!< The context material law (extracted from the spatialParams)

NEW_PROP_TAG(VtkAddDeltaP); // !< Output of pressure minus a fixed value
SET_BOOL_PROP(MPNC, VtkAddDeltaP, false);

//! The Model for kinetic Mass and Energy Transfer
NEW_PROP_TAG(MPNCModelKinetic);

//! Enable kinetic resolution of mass transfer processes?
NEW_PROP_TAG(EnableKinetic);

//! Property for the definition of the number of energy equations (0,1,2,3)
NEW_PROP_TAG(NumEnergyEquations);

//! average the velocity in the model
NEW_PROP_TAG(VelocityAveragingInModel);

//! which nusselt correlation to use
NEW_PROP_TAG(NusseltFormulation);

//! which sherwood correlation to use
NEW_PROP_TAG(SherwoodFormulation);

}
}

#endif
