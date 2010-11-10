// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Andreas Lauser                                    *
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
 * \brief Defines a type tags and some fundamental properties for
 *        fully coupled and decoupled models
 */
#ifndef DUMUX_BASIC_PROPERTIES_HH
#define DUMUX_BASIC_PROPERTIES_HH

#include <dumux/common/propertysystem.hh>

namespace Dumux
{
namespace Properties
{
///////////////////////////////////
// Type tag definitions:
//
// NumericModel
// |
// +-> CoupledModel
// |
// \-> DecoupledModel
///////////////////////////////////

//! Type tag for all models.
NEW_TYPE_TAG(NumericModel);

//! Type tag for all fully coupled models.
NEW_TYPE_TAG(CoupledModel, INHERITS_FROM(NumericModel));

//! Type tag for all decoupled models.
NEW_TYPE_TAG(DecoupledModel, INHERITS_FROM(NumericModel));


///////////////////////////////////
// Property names which are always available:
//
// Scalar
///////////////////////////////////

//! Property to specify the type of scalar values.
NEW_PROP_TAG(Scalar);

///////////////////////////////////
// Default values for properties:
//
// Scalar -> double
///////////////////////////////////

//! Set the default type of scalar values to double
SET_TYPE_PROP(NumericModel, Scalar, double);

} // namespace Properties
} // namespace Dumux

#endif
