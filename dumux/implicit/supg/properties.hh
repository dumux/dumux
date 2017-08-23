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
#ifndef DUMUX_FEM_PROPERTIES_HH
#define DUMUX_FEM_PROPERTIES_HH

#include <dumux/implicit/properties.hh>

/*!
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup BoxModel
 * \file
 * \brief Specify the shape functions, operator assemblers, etc
 *        used for the BoxModel.
 */
namespace Dumux
{

namespace Properties
{
// \{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the box-scheme
NEW_TYPE_TAG(FemModel, INHERITS_FROM(ImplicitBase));

NEW_PROP_TAG(FemIntegrationPointData); //!< Stores data on shape values and gradients at an integration point
NEW_PROP_TAG(SecondaryVariables); //!< Container to store primary/secondary variables
NEW_PROP_TAG(FeBasis); //!< The global finite element basis type
NEW_PROP_TAG(FemBasisOrder); //! The order of the finite element basis
NEW_PROP_TAG(FemQuadratureOrder); //! The order of the gaussian integral rule
}
}

// \}

#include "propertydefaults.hh"

#endif
