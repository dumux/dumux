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
 * \ingroup Properties
 * \file
 *
 * \brief Declares properties required for finite-volume models models.
 */

#ifndef DUMUX_FV_PROPERTIES_HH
#define DUMUX_FV_PROPERTIES_HH

#include <dune/istl/bvector.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/properties/grid.hh>
#include <dumux/common/properties/numericmodel.hh>

#include <dumux/implicit/gridvariables.hh>

namespace Dumux
{
namespace Properties
{
//! Type tag for finite-volume schemes.
NEW_TYPE_TAG(FiniteVolumeModel, INHERITS_FROM(GridProperties, NumericModel));

//! The grid variables
SET_TYPE_PROP(FiniteVolumeModel, GridVariables, GridVariables<TypeTag>);

//! The type of a solution for a whole element
SET_TYPE_PROP(FiniteVolumeModel, ElementSolutionVector, Dune::BlockVector<typename GET_PROP_TYPE(TypeTag, PrimaryVariables)>);

//! We do not store the FVGeometry by default
SET_BOOL_PROP(FiniteVolumeModel, EnableFVGridGeometryCache, false);

//! We do not store the volume variables by default
SET_BOOL_PROP(FiniteVolumeModel, EnableGlobalVolumeVariablesCache, false);

//! disable flux variables data caching by default
SET_BOOL_PROP(FiniteVolumeModel, EnableGlobalFluxVariablesCache, false);

} // namespace Properties
} // namespace Dumux

 #endif
