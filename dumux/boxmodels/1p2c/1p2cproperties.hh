// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Defines the properties required for the single-phase,
 *        two-component BOX model.
 */

#ifndef DUMUX_1P2C_PROPERTIES_HH
#define DUMUX_1P2C_PROPERTIES_HH

#include<dumux/boxmodels/common/boxproperties.hh>

namespace Dumux
{
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the isothermal single-phase, two-component problems
NEW_TYPE_TAG(BoxOnePTwoC, INHERITS_FROM(BoxModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(NumPhases);   //!< Number of fluid phases in the system
NEW_PROP_TAG(NumComponents);   //!< Number of fluid components in the system
NEW_PROP_TAG(OnePTwoCIndices); //!< Enumerations for the 1p2c models
NEW_PROP_TAG(SpatialParameters); //!< The type of the spatial parameters
NEW_PROP_TAG(FluidSystem); //!< Type of the multi-component relations
NEW_PROP_TAG(UpwindAlpha);   //!< The default value of the upwind parameter
NEW_PROP_TAG(EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG(PhaseIndex); //!< The index of the phase in the fluid system
NEW_PROP_TAG(Comp1Index); //!< The index of the first component in the fluid system
NEW_PROP_TAG(Comp2Index); //!< The index of the second component in the fluid system
}

}

#endif

