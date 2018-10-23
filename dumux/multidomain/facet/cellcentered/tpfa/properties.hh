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
 * \ingroup MultiDomain
 * \ingroup FacetCoupling
 * \brief Properties (and default properties) for all models using cell-centered
 *        finite volume scheme with TPFA together with coupling across the grid element facets
 * \note If n is the dimension of the lowest grid to be considered in the hierarchy,
 *       all problem type tags for the grids with the dimension m > n must inherit
 *       from these or other facet coupling properties (e.g. BoxFacetCouplingModel).
 */

#ifndef DUMUX_FACETCOUPLING_CC_TPFA_PROPERTIES_HH
#define DUMUX_FACETCOUPLING_CC_TPFA_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/multidomain/facet/cellcentered/upwindscheme.hh>
#include <dumux/multidomain/facet/cellcentered/localresidual.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/darcyslaw.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/fickslaw.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/fourierslaw.hh>

#include <dumux/porousmediumflow/fluxvariables.hh>

namespace Dumux {

namespace Properties {

//! Type tag for the cell-centered tpfa scheme with coupling to
//! another sub-domain living on the grid facets.
NEW_TYPE_TAG(CCTpfaFacetCouplingModel, INHERITS_FROM(CCTpfaModel));

//! Use the tpfa facet coupling-specific Darcy's law
SET_TYPE_PROP(CCTpfaFacetCouplingModel,
              AdvectionType,
              CCTpfaFacetCouplingDarcysLaw< typename GET_PROP_TYPE(TypeTag, Scalar),
                                            typename GET_PROP_TYPE(TypeTag, FVGridGeometry) >);

//! Use the tpfa facet coupling-specific Fick's law
SET_TYPE_PROP(CCTpfaFacetCouplingModel, MolecularDiffusionType, CCTpfaFacetCouplingFicksLaw<TypeTag>);

//! Use the tpfa facet coupling-specific Fourier's law
SET_TYPE_PROP(CCTpfaFacetCouplingModel, HeatConductionType, CCTpfaFacetCouplingFouriersLaw<TypeTag>);

//! Use the cc local residual for models with facet coupling
SET_TYPE_PROP(CCTpfaFacetCouplingModel, BaseLocalResidual, CCFacetCouplingLocalResidual<TypeTag>);

//! Per default, use the porous medium flow flux variables with the modified upwind scheme
SET_TYPE_PROP(CCTpfaFacetCouplingModel,
              FluxVariables,
              PorousMediumFluxVariables<TypeTag, CCFacetCouplingUpwindScheme<typename GET_PROP_TYPE(TypeTag, FVGridGeometry)>>);

} // namespace Properties
} // namespace Dumux

#endif
