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
 * \param dx the x-distance between the left and right centers
 * \param dy the y-distance between the left and right centers
 *
 */
template<class Scalar, class GlobalPosition, class FVElementGeometry>
std::array<Scalar,3> diffusionFlux(const Scalar waterDepthLeft,
                                   const Scalar waterDepthRight,
                                   Scalar velocityXLeft,
                                   Scalar velocityXRight,
                                   Scalar velocityYLeft,
                                   Scalar velocityYRight,
                                   const GlobalPosition& nxy,
                                   const Scalar dx,
                                   const Scalar dy,
                                   const FVElementGeometry& fvGeometry,
                                   const typename FVElementGeometry::SubControlVolumeFace& scvf)
{
    using std::abs;
    // Initialization
    Scalar gradU = 0.0;
    Scalar gradV = 0.0;

    Scalar distance = std::sqrt(dx*dx + dy*dy);

    // Compute the velocity gradients
    // Only if both cells have sufficient water
    // Does this result in problems with differentiability?
    if ((waterDepthLeft > 0.01) & (waterDepthRight > 0.01))
    {
        gradU = (velocityXRight-velocityXLeft)/distance;
        gradV = (velocityYRight-velocityYLeft)/distance;
    }

    // Now get the (constant) background turbulent viscosity
    static const Scalar turbBGViscosity = getParam<Scalar>("Problem.TurbViscosity", 1.0e-6);

    Scalar turbViscosityLeft  = turbBGViscosity;
    Scalar turbViscosityRight = turbBGViscosity;

    // Use a harmonic average of the viscosity and the depth at the interface.
    Scalar viscDepthAverage = 2.0*(turbViscosityLeft*waterDepthLeft*turbViscosityRight*waterDepthRight)/(turbViscosityLeft*waterDepthLeft + turbViscosityRight*waterDepthRight);

    // Factor that takes the direction of the unit vector into account
    Scalar fac = 2.0*(nxy[0]*dx + nxy[1]*dy)/distance;

    // Compute the diffusive fluxes
    Scalar uPartOfUDiffusion = viscDepthAverage * gradU * fac;
    Scalar vPartOfVDiffusion = viscDepthAverage * gradV * fac;

    // The actual stress term / diffusion flux for this edge
    Scalar uDiffusion = uPartOfUDiffusion;
    Scalar vDiffusion = vPartOfVDiffusion;

    std::array<Scalar, 3> localFlux;
    localFlux[0] = 0.0;
    localFlux[1] = -uDiffusion;
    localFlux[2] = -vDiffusion;

    return localFlux;
}

} // end namespace ShallowWater
} // end namespace Dumux

#endif
