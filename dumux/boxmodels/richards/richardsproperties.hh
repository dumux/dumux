// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Institute of Hydraulic Engineering                                      *
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
 * \file
 *
 * \brief Contains the property declarations for the Richards box
 *        model.
 */
#ifndef DUMUX_RICHARDS_PROPERTIES_HH
#define DUMUX_RICHARDS_PROPERTIES_HH

#include <dumux/boxmodels/common/boxproperties.hh>

namespace Dumux
{
/*!
 * \addtogroup RichardsModel
 */
// \{
///////////////////////////////////////////////////////////////////////////
// properties for the isothermal richards model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for problems discretized using the Richards model
NEW_TYPE_TAG(BoxRichards, INHERITS_FROM(BoxModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(RichardsIndices); //!< Enumerations used by the Richards models
NEW_PROP_TAG(SpatialParameters); //!< The type of the spatial parameters object
NEW_PROP_TAG(MaterialLaw);   //!< The material law which ought to be used (by default extracted from the spatial parameters)
NEW_PROP_TAG(MaterialLawParams); //!< The type of the parameter object for the material law (by default extracted from the spatial parameters)
NEW_PROP_TAG(FluidSystem); //!< The fluid system to be used for the Richards model
NEW_PROP_TAG(FluidState); //!< The fluid state to be used for the Richards model
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(MobilityUpwindAlpha); //!< The value of the upwind parameter for the mobility
// \}
};

} // end namepace

#endif
