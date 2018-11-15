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
 * \brief Properties (and default properties) for all models using the box
 *        scheme together with coupling across the grid element facets
 * \note If n is the dimension of the lowest grid to be considered in the hierarchy,
 *       all problem type tags for the grids with the dimension m > n must inherit
 *       from these or other facet coupling properties (e.g. CCTpfaFacetCouplingModel).
 */

#ifndef DUMUX_FACETCOUPLING_BOX_PROPERTIES_HH
#define DUMUX_FACETCOUPLING_BOX_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/box/properties.hh>

#include <dumux/multidomain/facet/box/darcyslaw.hh>
#include <dumux/multidomain/facet/box/elementboundarytypes.hh>
#include <dumux/multidomain/facet/box/fvgridgeometry.hh>
#include <dumux/multidomain/facet/box/localresidual.hh>
#include <dumux/multidomain/facet/box/upwindscheme.hh>

#include <dumux/porousmediumflow/fluxvariables.hh>

namespace Dumux {

namespace Properties {

//! Type tag for the box scheme with coupling to
//! another sub-domain living on the grid facets.
// Create new type tags
namespace TTag {
struct BoxFacetCouplingModel { using InheritsFrom = std::tuple<BoxModel>; };
} // end namespace TTag

//! Use the box local residual for models with facet coupling
SET_TYPE_PROP(BoxFacetCouplingModel, BaseLocalResidual, BoxFacetCouplingLocalResidual<TypeTag>);

//! Use the box facet coupling-specific Darcy's law
SET_TYPE_PROP(BoxFacetCouplingModel,
              AdvectionType,
              BoxFacetCouplingDarcysLaw< GetPropType<TypeTag, Properties::Scalar>,
                                         GetPropType<TypeTag, Properties::FVGridGeometry> >);

//! Per default, use the porous medium flow flux variables with the modified upwind scheme
SET_TYPE_PROP(BoxFacetCouplingModel,
              FluxVariables,
              PorousMediumFluxVariables<TypeTag, BoxFacetCouplingUpwindScheme<GetPropType<TypeTag, Properties::FVGridGeometry>>>);

//! Per default, use the porous medium flow flux variables with the modified upwind scheme
SET_TYPE_PROP(BoxFacetCouplingModel,
              ElementBoundaryTypes,
              BoxFacetCouplingElementBoundaryTypes<GetPropType<TypeTag, Properties::BoundaryTypes>>);

//! Set the default for the grid finite volume geometry
SET_PROP(BoxFacetCouplingModel, FVGridGeometry)
{
private:
    static constexpr bool enableCache = GET_PROP_VALUE(TypeTag, EnableFVGridGeometryCache);
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = BoxFacetCouplingFVGridGeometry<Scalar, GridView, enableCache>;
};

} // namespace Properties
} // namespace Dumux

#endif
