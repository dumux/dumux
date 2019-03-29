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
 *
 * \brief Defines the properties required for the richards linear-elastic model.
 *
 * This class inherits from the properties of the richards model and
 * from the properties of the simple linear-elastic model
 */

#ifndef DUMUX_ELASTICRICHARDS_PROPERTIES_HH
#define DUMUX_ELASTICRICHARDS_PROPERTIES_HH

#include <dumux/implicit/box/properties.hh>
#include <dumux/porousmediumflow/richards/implicit/properties.hh>

namespace Dumux
{
////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{
//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the richards model with a linear elastic matrix
NEW_TYPE_TAG(BoxElasticRichards, INHERITS_FROM(BoxModel));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG(DisplacementGridFunctionSpace); //!< grid function space for the displacement
NEW_PROP_TAG(PressureGridFunctionSpace); //!< grid function space for the pressure, saturation, ...
NEW_PROP_TAG(GridOperatorSpace); //!< The grid operator space
NEW_PROP_TAG(GridOperator); //!< The grid operator space
NEW_PROP_TAG(PressureFEM); //!< FE space used for pressure, saturation, ...
NEW_PROP_TAG(DisplacementFEM); //!< FE space used for displacement

//! Returns whether the output should be written according to
//! rock mechanics sign convention (compressive stresses > 0)
NEW_PROP_TAG(VtkRockMechanicsSignConvention);

//! Specifies the grid function space used for sub-problems
NEW_PROP_TAG(GridFunctionSpace);

//! Specifies the grid operator used for sub-problems
NEW_PROP_TAG(GridOperator);

//! Specifies the grid operator space used for sub-problems
NEW_PROP_TAG(GridOperatorSpace);

//! Specifies the type of the constraints
NEW_PROP_TAG(Constraints);

//! Specifies the type of the constraints transformation
NEW_PROP_TAG(ConstraintsTrafo);

//! Specifies the local finite element space
NEW_PROP_TAG(LocalFEMSpace);

//! Specifies the local operator
NEW_PROP_TAG(LocalOperator);

//! The type traits required for using the AMG backend
NEW_PROP_TAG(AmgTraits);

NEW_PROP_TAG(EffectivePermeabilityModel);

NEW_PROP_TAG(UseHead); //!< Defines whether pressure [Pa] (false) or pressure head [cm] (true) is used. got from richards properties
}

}

#endif
