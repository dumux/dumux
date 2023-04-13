// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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
 * \brief Prepare and compute the shallow water advective flux.
 *
 * Prepares the Riemann problem for the advective flux for the 2D shallow
 * water model. The actual model uses an exact Riemann solver after Torro
 * and the reconstruction after Audusse. A flux limiter is
 * applied to limit water flow for small water depths.
 *
 * The computed water flux of the Riemann solver is given in m^2/s, the
 * momentum fluxes are given in m^3/s^2. The Riemann flux is multiplied by
 * scvf.area() (given in m) to obtain the flux over the face.
 *
 * \todo Add more numerical fluxes and reconstruction methods.
 */
template<class NumEqVector>
class ShallowWaterFlux
{

public:

    using Cache = FluxVariablesCaching::EmptyAdvectionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller;

    /*!
     * \ingroup Flux
     * \brief Prepares and compute the shallow water advective flux.
     *
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
