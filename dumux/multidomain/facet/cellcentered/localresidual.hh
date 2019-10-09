// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup FacetCoupling
 * \brief Calculates the element-wise residual for cell-centered discretization schemes
 *        in models where coupling occurs across the element facets. This extra implementation
 *        is necessary as facets that lie on the boundary but couple to a facet element have to be
 *        treated differently.
 */
#ifndef DUMUX_FACETCOUPLING_CC_LOCAL_RESIDUAL_HH
#define DUMUX_FACETCOUPLING_CC_LOCAL_RESIDUAL_HH

#include <dumux/common/properties.hh>
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

    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;

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
