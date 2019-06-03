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
 * \brief This file includes a function to to compute the diffusive flux
 *
 */
#ifndef DUMUX_FLUX_SHALLOW_WATER_DIFFUSION_FLUX_HH
#define DUMUX_FLUX_SHALLOW_WATER_DIFFUSION_FLUX_HH

#include <dumux/flux/shallowwater/fluxlimiterlet.hh>

namespace Dumux {
namespace ShallowWater {

/*!
 * \ingroup ShallowWaterFlux
 * \brief Compute the diffusive flux contribution from the interface
 *        shear stress
 *
 *
 * The diffusive flux
 * \f[
 * \int \int_{V} \mathbf{\nabla} \cdot \nu_t h \mathbf{\nabla} \mathbf{u} dV
 * \f]
 * is re-written using Gauss' divergence theorem to:
 * \f[
 * \int_{S_f} \nu_t h \mathbf{\nabla} \mathbf{u} \cdot \mathbf{n_f} dS
 * \f]
 *
 * \param waterDepthLeft water depth on the left side
 * \param waterDepthRight water depth on the right side
 * \param velocityXLeft veloctiyX on the left side
 * \param velocityXRight velocityX on the right side
 * \param velocityYLeft velocityY on the left side
 * \param velocityYRight velocityY on the right side
 * \param nxy the normal vector
 * \param distance the distance between the left and right centers
 *
 */
template<class Scalar, class GlobalPosition>
std::array<Scalar,3> diffusionFlux(const Scalar waterDepthLeft,
                                   const Scalar waterDepthRight,
                                   Scalar velocityXLeft,
                                   Scalar velocityXRight,
                                   Scalar velocityYLeft,
                                   Scalar velocityYRight,
                                   const GlobalPosition& nxy,
                                   const Scalar distance)
{
    using std::max;

    // Initialization
    Scalar uPartOfUDiffusion = 0.0;
    Scalar vPartOfUDiffusion = 0.0;
    Scalar uPartOfVDiffusion = 0.0;
    Scalar vPartOfVDiffusion = 0.0;
    Scalar gradU = 0.0;
    Scalar gradV = 0.0;

    // TODO: Get the viscosity from a specified turbulence model
    // e.g. getTurbulentViscosity()
    // For now use the kinematic viscosity \nu = 1.0e-6
    Scalar kinematicViscosity = 1.0e-6;

    // The left and right (turbulent eddy) viscosities
    Scalar viscosityL = kinematicViscosity;
    Scalar viscosityR = kinematicViscosity;

    // Compute the velocity gradients
    // Only if both cells have sufficient water
    // Does this result in problems with differentiability?
    if ((waterDepthLeft > 0.01) & (waterDepthRight > 0.01)){
        gradU = (velocityXRight-velocityXLeft)/distance;
        gradV = (velocityYRight-velocityYLeft)/distance;
    }

    // For now use an arithmetic average of the viscosity and the depth at the interface.
    // Perhaps a harmonic average is better suited for strongly varying viscosity or depth
    Scalar viscDepthAverage = 0.5*(viscosityL*waterDepthLeft + viscosityR*waterDepthRight);

    // Compute the diffusive fluxes
    uPartOfUDiffusion = viscDepthAverage * gradU * (1.0 + nxy[0]*nxy[0]);
    vPartOfUDiffusion = viscDepthAverage * gradV * nxy[0]*nxy[1];
    uPartOfVDiffusion = viscDepthAverage * gradU * nxy[0]*nxy[1];
    vPartOfVDiffusion = viscDepthAverage * gradV * (1.0 + nxy[1]*nxy[1]);


    // The actual stress term / diffusion flux for this edge
    Scalar uDiffusion = uPartOfUDiffusion + vPartOfUDiffusion;
    Scalar vDiffusion = uPartOfVDiffusion + vPartOfVDiffusion;

    // Compute the mobility of the flux with the fluxlimiter
    /*const Scalar mobility = ShallowWater::fluxLimiterLET(waterDepthLeft, //Reconstructed,
                                                         waterDepthRight, //Reconstructed,
                                                         0.001,
                                                         0.00001);
    */

    std::array<Scalar, 3> localFlux;
    localFlux[0] = 0.0;
    localFlux[1] = -uDiffusion; // * mobility;
    localFlux[2] = -vDiffusion; // * mobility;

    return localFlux;
}

} // end namespace ShallowWater
} // end namespace Dumux

#endif
