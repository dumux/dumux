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
 * \ingroup MixedDimension
 * \brief Base properties for the bulk problems in mixed dimensional models
 *        with a lower dimensional model living on the element facets.
 */

#ifndef DUMUX_FACET_MIXEDDIMENSION_PROPERTIES_HH
#define DUMUX_FACET_MIXEDDIMENSION_PROPERTIES_HH

#include <dumux/mixeddimension/properties.hh>
#include <dumux/mixeddimension/facet/subproblemlocaljacobian.hh>

namespace Dumux
{

namespace Properties
{
//! The type tag for models using facet coupling (only overwrites the local jacobians)
NEW_TYPE_TAG(MixedDimensionFacetCoupling, INHERITS_FROM(MixedDimension));
SET_TYPE_PROP(MixedDimensionFacetCoupling, BulkLocalJacobian, FacetCouplingBulkLocalJacobian<TypeTag>);
SET_TYPE_PROP(MixedDimensionFacetCoupling, LowDimLocalJacobian, FacetCouplingLowDimLocalJacobian<TypeTag>);

}//end namespace Properties

}//end namespace Dumux

#endif
