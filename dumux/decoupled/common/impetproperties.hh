// $Id$
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
#ifndef DUMUX_IMPET_PROPERTIES_HH
#define DUMUX_IMPET_PROPERTIES_HH

#include <dumux/decoupled/common/decoupledproperties.hh>
/*!
 * \ingroup Sequential
 * \defgroup IMPET IMPET (Implicit Pressure Explicit Transport)
 */
/*!
 * \ingroup IMPET
 * \ingroup Properties
 */
/*!
 * \file
 * \brief Base file for properties related to sequential IMPET algorithms
 */
namespace Dumux
{

template<class TypeTag>
class IMPET;

namespace Properties
{
/*!
 *
 * \brief General properties for sequential IMPET algorithms
 *
 * This class holds properties necessary for the sequential IMPET solution.
 */

//////////////////////////////////////////////////////////////////
// Type tags tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the diffusion-scheme
NEW_TYPE_TAG(IMPET, INHERITS_FROM(DecoupledModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG(PressureModel);         //!< The type of the discretizations
NEW_PROP_TAG(TransportModel);         //!< The type of the discretizations

NEW_PROP_TAG(CFLFactor);
NEW_PROP_TAG(IterationFlag);
NEW_PROP_TAG(IterationNumber);
NEW_PROP_TAG(MaximumDefect);
NEW_PROP_TAG(RelaxationFactor);

SET_TYPE_PROP(IMPET, Model, IMPET<TypeTag>);

SET_SCALAR_PROP(IMPET, CFLFactor, 1);
SET_INT_PROP(IMPET, IterationFlag, 0);
SET_INT_PROP(IMPET, IterationNumber, 2);
SET_SCALAR_PROP(IMPET, MaximumDefect, 1e-5);
SET_SCALAR_PROP(IMPET, RelaxationFactor, 1);

}
}

#endif
