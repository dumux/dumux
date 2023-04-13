// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
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
#include <dumux/discretization/cctpfa.hh>
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
// Create new type tags
namespace TTag {
struct CCTpfaFacetCouplingModel { using InheritsFrom = std::tuple<CCTpfaModel>; };
} // end namespace TTag

//! Use the tpfa facet coupling-specific Darcy's law
template<class TypeTag>
struct AdvectionType<TypeTag, TTag::CCTpfaFacetCouplingModel>
{
    using type = CCTpfaFacetCouplingDarcysLaw< GetPropType<TypeTag, Properties::Scalar>,
                                               GetPropType<TypeTag, Properties::GridGeometry> >;
};

//! Use the tpfa facet coupling-specific Ficks's law
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::CCTpfaFacetCouplingModel>
{
    using type = CCTpfaFacetCouplingFicksLaw< TypeTag >;
};

//! Use the tpfa facet coupling-specific Ficks's law
template<class TypeTag>
struct HeatConductionType<TypeTag, TTag::CCTpfaFacetCouplingModel>
{
    using type = CCTpfaFacetCouplingFouriersLaw< TypeTag >;
};

//! Use the cc local residual for models with facet coupling
template<class TypeTag>
struct BaseLocalResidual<TypeTag, TTag::CCTpfaFacetCouplingModel> { using type = CCFacetCouplingLocalResidual<TypeTag>; };

//! Per default, use the porous medium flow flux variables with the modified upwind scheme
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::CCTpfaFacetCouplingModel>
{
    using type = PorousMediumFluxVariables<TypeTag,
                                           CCFacetCouplingUpwindScheme<GetPropType<TypeTag, Properties::GridGeometry>>>;
};

} // namespace Properties
} // namespace Dumux

#endif
