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
 * \ingroup ShallowWaterFlux
 * \brief This file includes a function to construct the Riemann problem
 *
 */
#ifndef DUMUX_FLUX_SHALLOW_WATER_RIEMANN_PROBLEM_HH
#define DUMUX_FLUX_SHALLOW_WATER_RIEMANN_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/flux/shallowwater/fluxlimiterlet.hh>
#include <dumux/flux/shallowwater/exactriemann.hh>

namespace Dumux {
namespace ShallowWater {

/*!
 * \ingroup ShallowWaterFlux
 * \brief Construct a Riemann problem and solve it
 *
 *
 * Riemann problem applies the hydrostatic reconstruction, uses the
 * Riemann invariants to transform the two-dimensional problem to an
 * one-dimensional problem and solves this new problem, and rotates
 * the problem back. Further it applies an flux limiter for the water
 * flux handle drying of elements.
 * The correction of the bed slope source term leads to a
 * non-symmetric flux term at the interface for the momentum equations
 * since DuMuX computes the fluxes twice from each side this does not
 * matter.
 *
 * So far we have only the exact Riemann solver, and the reconstruction
 * after Audusse but further solvers and reconstructions ca be
 * implemented.
 *
 * The computed water flux (localFlux[0]) is given in m^2/s, the
 * momentum fluxes (localFlux[1], localFlux[2]) are given in m^3/s^2.
 * Later this flux will be multiplied by the scvf.area() (given in m
 * for a 2D problem) to get the flux over a face.
 *
 * \param waterDepthLeft water depth on the left side
 * \param waterDepthRight water depth on the right side
 * \param velocityXLeft veloctiyX on the left side
 * \param velocityXRight velocityX on the right side
 * \param velocityYLeft velocityY on the left side
 * \param velocityYRight velocityY on the right side
 * \param bedSurfaceLeft surface of the bed on the left side
 * \param bedSurfaceRight surface of the bed on the right side
 * \param gravity gravity constant
 * \param nxy the normal vector
 *
 */
template<class Scalar, class GlobalPosition>
std::array<Scalar,3> riemannProblem(const Scalar waterDepthLeft,
                                    const Scalar waterDepthRight,
                                    Scalar velocityXLeft,
                                    Scalar velocityXRight,
                                    Scalar velocityYLeft,
                                    Scalar velocityYRight,
                                    const Scalar bedSurfaceLeft,
                                    const Scalar bedSurfaceRight,
                                    const Scalar gravity,
                                    const GlobalPosition& nxy)
{
    using std::max;

    // hydrostatic reconstrucion after Audusse
    const Scalar dzl = max(0.0, bedSurfaceRight - bedSurfaceLeft);
    const Scalar waterDepthLeftReconstructed = max(0.0, waterDepthLeft - dzl);
    const Scalar dzr = max(0.0, bedSurfaceLeft - bedSurfaceRight);
    const Scalar waterDepthRightReconstructed = max(0.0, waterDepthRight - dzr);

    // make rotation of the flux we compute an 1d flux
    Scalar tempFlux = velocityXLeft;
    velocityXLeft =  nxy[0] * tempFlux + nxy[1] * velocityYLeft;
    velocityYLeft = -nxy[1] * tempFlux + nxy[0] * velocityYLeft;

    tempFlux = velocityXRight;
    velocityXRight =  nxy[0] * tempFlux + nxy[1] * velocityYRight;
    velocityYRight = -nxy[1] * tempFlux + nxy[0] * velocityYRight;

    auto riemannResult = ShallowWater::exactRiemann(waterDepthLeftReconstructed,
                                                    waterDepthRightReconstructed,
                                                    velocityXLeft,
                                                    velocityXRight,
                                                    velocityYLeft,
                                                    velocityYRight,
                                                    gravity);

    //redo rotation
    tempFlux = riemannResult.flux[1];
    riemannResult.flux[1] = nxy[0] * tempFlux - nxy[1] * riemannResult.flux[2];
    riemannResult.flux[2] = nxy[1] * tempFlux + nxy[0] * riemannResult.flux[2];

    // Add reconstruction flux from Audusse reconstruction
    const Scalar hgzl = 0.5 * (waterDepthLeftReconstructed + waterDepthLeft) * (waterDepthLeftReconstructed - waterDepthLeft);
    const Scalar hdxzl = gravity * nxy[0] * hgzl;
    const Scalar hdyzl = gravity * nxy[1] * hgzl;

    /*Right side is computed from the other side otherwise the
      following "non-symetric" fluxes are needed:

      Scalar hgzr = 0.5 * (waterDepthRightReconstructed + waterDepthRight) * (waterDepthRightReconstructed - waterDepthRight);
      Scalar hdxzr = gravity * nxy[0] * hgzr;
      Scalar hdyzrhdyzr = gravity * nxy[1] * hgzr;
    */

    // compute the mobility of the flux with the fluxlimiter
    static const Scalar upperWaterDepthFluxLimiting = getParam<Scalar>("FluxLimiterLET.UpperWaterDepth", 1e-3);
    static const Scalar lowerWaterDepthFluxLimiting = getParam<Scalar>("FluxLimiterLET.LowerWaterDepth", 1e-5);
    static const bool upwindWaterDepthFluxLimiting = getParam<bool>("FluxLimiterLET.UpwindFluxLimiting", false);

    Scalar limitingDepth = (waterDepthLeftReconstructed + waterDepthRightReconstructed) * 0.5;

    //Using the upwind water depth from the flux direction can improve stability.
    if (upwindWaterDepthFluxLimiting)
    {
        if (riemannResult.flux[0] < 0)
        {
            limitingDepth = waterDepthRightReconstructed;
        }else
        {
            limitingDepth = waterDepthLeftReconstructed;
        }
    }

    const Scalar mobility = ShallowWater::fluxLimiterLET(limitingDepth,
                                                         limitingDepth,
                                                         upperWaterDepthFluxLimiting,
                                                         lowerWaterDepthFluxLimiting);
    std::array<Scalar, 3> localFlux;
    localFlux[0] = riemannResult.flux[0] * mobility;
    localFlux[1] = (riemannResult.flux[1] - hdxzl);
    localFlux[2] = (riemannResult.flux[2] - hdyzl);

    return localFlux;
}

} // end namespace ShallowWater
} // end namespace Dumux

#endif
