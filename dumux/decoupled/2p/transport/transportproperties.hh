// $Id: transportproperties.hh 3732 2010-06-11 13:27:20Z bernd $
/*****************************************************************************
 *   Copyright (C) 2009 by Markus Wolff                                      *
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
#ifndef DUMUX_TRANSPORT_PROPERTIES_HH
#define DUMUX_TRANSPORT_PROPERTIES_HH


/*!
 * \file
 * \brief Specify the shape functions, operator assemblers, etc
 *        used for the BoxScheme.
 */
namespace Dumux
{

template<class TypeTag>
class DiffusivePart;

template<class TypeTag>
class ConvectivePart;

namespace Properties
{
/*!
 * \addtogroup diffusion
 */
// \{

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(Transport);

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG( DiffusivePart );         //!< The type of the diffusive part in a transport equation
NEW_PROP_TAG( ConvectivePart );        //!< The type of a convective part in a transport equation

SET_TYPE_PROP(Transport, DiffusivePart, DiffusivePart<TypeTag>);
SET_TYPE_PROP(Transport, ConvectivePart, ConvectivePart<TypeTag>);

}
}

#endif
