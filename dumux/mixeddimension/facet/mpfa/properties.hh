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

#include <dumux/implicit/cellcentered/mpfa/properties.hh>

#include <dumux/mixeddimension/subproblemproperties.hh>
#include <dumux/mixeddimension/facet/mpfa/interactionvolume.hh>
#include <dumux/mixeddimension/facet/mpfa/interiorboundarydata.hh>
#include <dumux/mixeddimension/facet/mpfa/darcyslaw.hh>
#include <dumux/mixeddimension/facet/mpfa/fickslaw.hh>
#include <dumux/mixeddimension/facet/mpfa/fourierslaw.hh>

namespace Dumux
{

namespace Properties
{
NEW_TYPE_TAG(FacetCouplingBulkMpfaModel, INHERITS_FROM(CCMpfaModel));

//! The boundary interaction volume class (we use the facet coupling specialized o-method interaction volume)
SET_TYPE_PROP(FacetCouplingBulkMpfaModel, BoundaryInteractionVolume, CCMpfaOFacetCouplingInteractionVolume<TypeTag>);

//! The interior boundary data class
SET_TYPE_PROP(FacetCouplingBulkMpfaModel, InteriorBoundaryData, CCMpfaFacetCouplingInteriorBoundaryData<TypeTag>);

//! Ẃe always enable interior boundaries
SET_BOOL_PROP(FacetCouplingBulkMpfaModel, EnableInteriorBoundaries, true);

//! Facet coupling is always true here
SET_BOOL_PROP(FacetCouplingBulkMpfaModel, MpfaFacetCoupling, true);

//! Darcy's Law
SET_TYPE_PROP(FacetCouplingBulkMpfaModel, AdvectionType, CCMpfaFacetCouplingDarcysLaw<TypeTag>);

//! Ficks's Law
SET_TYPE_PROP(FacetCouplingBulkMpfaModel, MolecularDiffusionType, CCMpfaFacetCouplingFicksLaw<TypeTag>);

//! Fourier's Law
SET_TYPE_PROP(FacetCouplingBulkMpfaModel, HeatConductionType, CCMpfaFacetCouplingFouriersLaw<TypeTag>);

}//end namespace Properties

}//end namespace Dumux

#endif
