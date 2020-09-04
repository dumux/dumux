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
 * \copydoc Dumux::ShallowWaterFlux
 */
#ifndef DUMUX_FLUX_SHALLOW_WATER_FLUX_HH
#define DUMUX_FLUX_SHALLOW_WATER_FLUX_HH

#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/flux/shallowwater/riemannproblem.hh>

namespace Dumux {

/*!
 * \ingroup Flux
 * \brief Computes the shallow water flux by solving a riemann problem.
 */
template<class NumEqVector>
class ShallowWaterFlux
{

public:

    using Cache = FluxVariablesCaching::EmptyAdvectionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller;

    /*!
     * \ingroup Flux
     * \brief Prepares the Riemann problem for the advective flux for
     *        the 2D shallow water model. The actual model uses an
     *        exact Riemann solver after Torro and the reconstruction
     *        after Audusse and a flux limiter for small water depths.
     *
     *        The computed water flux of the Riemann solver is given
     *        in m^2/s, the momentum fluxes are given in m^3/s^2. The
     *        Riemann flux is multiplied by scvf.area() (given in m
     *        for a 2D domain) to get the flux over the face.
     *
     * \todo The choice of the Riemann solver should be more flexible
     */
    template<class Problem, class FVElementGeometry, class ElementVolumeVariables>
    static NumEqVector flux(const Problem& problem,
                            const typename FVElementGeometry::GridGeometry::GridView::template Codim<0>::Entity& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const typename FVElementGeometry::SubControlVolumeFace& scvf)
    {

        //Get the inside and outside volume variables
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto& nxy = scvf.unitOuterNormal();
        const auto gravity = problem.spatialParams().gravity(scvf.center());

        auto riemannFlux = ShallowWater::riemannProblem(insideVolVars.waterDepth(),
                                                        outsideVolVars.waterDepth(),
                                                        insideVolVars.velocity(0),
                                                        outsideVolVars.velocity(0),
                                                        insideVolVars.velocity(1),
                                                        outsideVolVars.velocity(1),
                                                        insideVolVars.bedSurface(),
                                                        outsideVolVars.bedSurface(),
                                                        gravity,
                                                        nxy);

        NumEqVector localFlux(0.0);
        localFlux[0] = riemannFlux[0] * scvf.area();
        localFlux[1] = riemannFlux[1] * scvf.area();
        localFlux[2] = riemannFlux[2] * scvf.area();

        return localFlux;
    }
};

} // end namespace Dumux

#endif
