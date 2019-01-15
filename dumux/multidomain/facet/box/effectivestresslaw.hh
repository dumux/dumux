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
 * \brief Specialization of the effective stress law for the box scheme in the
 *        context of coupling across the bulk domain elemt facets. This computes
 *        the stress tensor and surface forces resulting from mechanical deformation
 *        and the pore pressure.
 */
#ifndef DUMUX_FACETCOUPLING_BOX_EFFECTIVE_STRESS_LAW_HH
#define DUMUX_FACETCOUPLING_BOX_EFFECTIVE_STRESS_LAW_HH

#include <dumux/flux/box/effectivestresslaw.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup FacetCoupling
 * \brief Effective stress law for the box facet coupling scheme
 * \tparam StressType type used for the computation of
 *         purely mechanical stresses (i.e. material law)
 * \tparam FVGridGeometry the finite volume grid geometry
 */
template<class StressType, class FVGridGeometry>
class BoxFacetCouplingEffectiveStressLaw
: public EffectiveStressLaw<StressType, FVGridGeometry, DiscretizationMethod::box>
{
    using ParentType = EffectiveStressLaw<StressType, FVGridGeometry, DiscretizationMethod::box>;

    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    //! export the type used for force vectors
    using typename ParentType::ForceVector;

    //! computes the force acting on a sub-control volume face
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static ForceVector force(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace& scvf,
                             const ElementFluxVarsCache& elemFluxVarCache)
    {
        if (!scvf.interiorBoundary())
            return ParentType::force(problem, element, fvGeometry, elemVolVars, scvf, elemFluxVarCache);
        else
            return problem.interiorBoundaryForce(element, fvGeometry, elemVolVars, scvf);
    }
};

} // end namespace Dumux

#endif // DUMUX_FACETCOUPLING_BOX_EFFECTIVE_STRESS_LAW_HH
