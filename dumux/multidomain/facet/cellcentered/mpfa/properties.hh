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
 *        finite volume scheme with MPFA together with coupling across the grid element facets
 * \note If n is the dimension of the lowest grid to be considered in the hierarchy,
 *       all problem type tags for the grids with the dimension m > n must inherit
 *       from these or other facet coupling properties (e.g. BoxFacetCouplingModel).
 */

#ifndef DUMUX_FACETCOUPLING_CC_MPFA_PROPERTIES_HH
#define DUMUX_FACETCOUPLING_CC_MPFA_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/ccmpfa.hh>

#include <dumux/multidomain/facet/cellcentered/upwindscheme.hh>
#include <dumux/multidomain/facet/cellcentered/localresidual.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/interactionvolume.hh>

#include <dumux/porousmediumflow/fluxvariables.hh>

namespace Dumux {

namespace Properties {

//! Type tag for the cell-centered mpfa scheme with coupling to
//! another sub-domain living on the grid facets.
// Create new type tags
namespace TTag {
struct CCMpfaFacetCouplingModel { using InheritsFrom = std::tuple<CCMpfaModel>; };
} // end namespace TTag

//! Use the cc local residual for models with facet coupling
template<class TypeTag>
struct BaseLocalResidual<TypeTag, TTag::CCMpfaFacetCouplingModel> { using type = CCFacetCouplingLocalResidual<TypeTag>; };

//! Use the facet coupling-specific mpfa-o interaction volume
template<class TypeTag>
struct PrimaryInteractionVolume<TypeTag, TTag::CCMpfaFacetCouplingModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NodalIndexSet = GetPropType<TypeTag, Properties::DualGridNodalIndexSet>;

    // use the default traits
    using Traits = CCMpfaOFacetCouplingDefaultInteractionVolumeTraits< NodalIndexSet, Scalar >;
public:
    using type = CCMpfaOFacetCouplingInteractionVolume< Traits >;
};

//! Use the facet coupling-specific mpfa-o interaction volume
template<class TypeTag>
struct SecondaryInteractionVolume<TypeTag, TTag::CCMpfaFacetCouplingModel>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NodalIndexSet = GetPropType<TypeTag, Properties::DualGridNodalIndexSet>;

    // use the default traits
    using Traits = CCMpfaOFacetCouplingDefaultInteractionVolumeTraits< NodalIndexSet, Scalar >;
public:
    using type = CCMpfaOFacetCouplingInteractionVolume< Traits >;
};

//! Per default, use the porous medium flow flux variables with the modified upwind scheme
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::CCMpfaFacetCouplingModel>
{
    using type = PorousMediumFluxVariables<TypeTag,
                                           CCFacetCouplingUpwindScheme<GetPropType<TypeTag, Properties::GridGeometry>>>;
};

} // namespace Properties
} // namespace Dumux

#endif
