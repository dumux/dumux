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
#ifndef DUMUX_CC_MPFA_PROPERTIES_HH
#define DUMUX_CC_MPFA_PROPERTIES_HH

#include <dumux/implicit/properties.hh>
#include <dumux/implicit/cellcentered/properties.hh>

/*!
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup CCMpfaModel
 * \file
 * \brief Specify classes used in Mpfa models.
 */
namespace Dumux
{

namespace Properties
{
// \{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for models based on the cell-centered multi-point flux approximation scheme
NEW_TYPE_TAG(CCMpfaModel, INHERITS_FROM(CCModel));

NEW_PROP_TAG(MpfaMethod); //! Specifies the mpfa method to be used
NEW_PROP_TAG(MpfaHelper); //! A Helper class depending on the mpfa method and grid dimension
NEW_PROP_TAG(InteractionVolume); //! The inner interaction volume type
NEW_PROP_TAG(BoundaryInteractionVolume); //! The interaction volume type used on the boundaries
NEW_PROP_TAG(GlobalInteractionVolumeSeeds); //! Class storing and managing the interaction volume seeds
NEW_PROP_TAG(QuadraturePoint); //! The quadrature point on the sub control volume faces (0.0 <= q <= 1.0)
NEW_PROP_TAG(FacetCoupling); //! This property has to be set to true if a sub model on the element facets should be existing
NEW_PROP_TAG(UseTpfaBoundary); //! This property specifies whether or not tpfa is to be used to handle the boundary fluxes
NEW_PROP_TAG(Xi); //! Parameter used for interior neumann boundaries or coupling terms on the element facets (0.0 <= xi <= 1.0)
}
}

// \}

#include <dumux/implicit/cellcentered/mpfa/propertydefaults.hh>

#endif
