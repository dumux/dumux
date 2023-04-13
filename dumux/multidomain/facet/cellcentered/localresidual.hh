// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetCoupling
 * \brief Calculates the element-wise residual for cell-centered discretization schemes
 *        in models where coupling occurs across the element facets. This extra implementation
 *        is necessary as facets that lie on the boundary but couple to a facet element have to be
 *        treated differently.
 */
#ifndef DUMUX_FACETCOUPLING_CC_LOCAL_RESIDUAL_HH
#define DUMUX_FACETCOUPLING_CC_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/assembly/cclocalresidual.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Calculates the element-wise residual for cell-centered discretization schemes
 *        in models where coupling occurs across the element facets. We only overwrite the
 *        function for the computation of a flux across a single sub-control volume face,
 *        as we need to additionally check if a boundary face couples to a facet element.
 */
template<class TypeTag>
class CCFacetCouplingLocalResidual : public CCLocalResidual<TypeTag>
{
    using ParentType = CCLocalResidual<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

public:
    //! pull up the parent's constructor
    using ParentType::ParentType;
    //! export the type used for element residuals
    using ElementResidualVector = typename ParentType::ElementResidualVector;

    //! evaluate the flux residual for a sub control volume face
    using ParentType::evalFlux;
    template< class Problem >
    NumEqVector evalFlux(const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const ElementFluxVariablesCache& elemFluxVarsCache,
                         const SubControlVolumeFace& scvf) const
    {
        // Even if scvf.boundary=true, compute flux on interior boundaries
        if (problem.couplingManager().isOnInteriorBoundary(element, scvf))
            return this->asImp().computeFlux(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
        else
            return ParentType::evalFlux(problem, element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
    }
};

} // end namespace Dumux

#endif
