// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Klaus Mosthaf                                     *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008 Bernd Flemisch                                       *
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
 * \brief Defines the properties required for the non-isotherm two-phase BOX model.
 */
#ifndef DUMUX_2PNI_PROPERTIES_HH
#define DUMUX_2PNI_PROPERTIES_HH

#include <dumux/boxmodels/2p/2pproperties.hh>

namespace Dumux
{
/*!
 * \addtogroup TwoPNIBoxModel
 */
// \{

namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the non-isothermal two-phase problems
NEW_TYPE_TAG(BoxTwoPNI, INHERITS_FROM(BoxTwoP));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(TwoPNIIndices); //!< Enumerations for the non-isothermal 2p models
}
}

#endif
