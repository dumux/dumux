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
 * \ingroup Flux
 * \copydoc Dumux::ShallowWaterDiffusiveFlux
 */
#ifndef DUMUX_FLUX_SHALLOW_WATER_DIFFUSIVE_FLUX_HH
#define DUMUX_FLUX_SHALLOW_WATER_DIFFUSIVE_FLUX_HH

#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/flux/shallowwater/diffusionflux.hh>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Computes the shallow water diffusive flux by adding all surrounding shear stresses.
 */
template<class NumEqVector>
class ShallowWaterDiffusiveFlux
{

public:

    using Cache = FluxVariablesCaching::EmptyDiffusionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller;

    /*!
     * \ingroup Flux
     * \brief Prepares the riemann problem for the diffusive flux for
     *        the shallow water model. The actual model adds the
     *        contributionof the shear stress tensor between two
     *        control volumes and a mechanism to avoid degenerate
     *        expressions for small water depths.
     *
     * \todo The implementation now is the simplest one involving
     *       only direct neighbours. A more complete implementation
     *       with a more elaborate stencil that also takes into account
     *       the non-orthogonal contributions can be considered.
     */
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables>
    static NumEqVector flux(const Problem& problem,
                            const typename FVElementGeometry::FVGridGeometry::GridView::template Codim<0>::Entity& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {

        //Get the inside and outside volume variables
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto& nxy = scvf.unitOuterNormal();

        // Inside and outside subcontrolvolumes
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());

        // Distance between the two cell centers
        const auto distance = (outsideScv.center() - insideScv.center()).two_norm();

        auto diffusiveFlux = ShallowWater::diffusionFlux(insideVolVars.waterDepth(),
                                                         outsideVolVars.waterDepth(),
                                                         insideVolVars.velocity(0),
                                                         outsideVolVars.velocity(0),
                                                         insideVolVars.velocity(1),
                                                         outsideVolVars.velocity(1),
                                                         nxy,
                                                         distance);

        NumEqVector localFlux(0.0);
        localFlux[0] = diffusiveFlux[0] * scvf.area();
        localFlux[1] = diffusiveFlux[1] * scvf.area();
        localFlux[2] = diffusiveFlux[2] * scvf.area();

        return localFlux;
    }
};

} // end namespace Dumux

#endif
