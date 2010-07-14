// $Id: impesproperties.hh 3784 2010-06-24 13:43:57Z bernd $
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
#ifndef DUMUX_IMPES_PROPERTIES_HH
#define DUMUX_IMPES_PROPERTIES_HH

#include <dumux/decoupled/common/decoupledproperties.hh>

/*!
 * \file
 * \brief Specify the shape functions, operator assemblers, etc
 *        used for the BoxScheme.
 */
namespace Dumux
{

template<class TypeTag>
class IMPES;

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
NEW_TYPE_TAG(IMPES, INHERITS_FROM(DecoupledModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(PressureModel);         //!< The type of the discretizations
NEW_PROP_TAG(SaturationModel);         //!< The type of the discretizations

NEW_PROP_TAG(CFLFactor);
NEW_PROP_TAG(IterationFlag);
NEW_PROP_TAG(IterationNumber);
NEW_PROP_TAG(MaximumDefect);
NEW_PROP_TAG(RelaxationFactor);

SET_TYPE_PROP(IMPES, Model, IMPES<TypeTag>);

SET_SCALAR_PROP(IMPES, CFLFactor, 1);
SET_INT_PROP(IMPES, IterationFlag, 0);
SET_INT_PROP(IMPES, IterationNumber, 2);
SET_SCALAR_PROP(IMPES, MaximumDefect, 1e-5);
SET_SCALAR_PROP(IMPES, RelaxationFactor, 1);

}
}

#endif
