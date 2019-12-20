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
 * \ingroup ShallowWaterModel
 * \brief Compute boundary conditions for the Riemann Solver
 *
 * The boundary conditions are given at the the outer face of the
 * the boundary cells. In this form the boundary condition can't be
 * processed by the riemann Solver, because it needs two cell states, one at
 * each side of a face. Therefore the Riemann invariants are used to
 * calculate a virtual outer state.
 */
#ifndef DUMUX_SHALLOWWATER_BOUNDARYFLUXES_HH
#define DUMUX_SHALLOWWATER_BOUNDARYFLUXES_HH

#include <array>
#include <cmath>

namespace Dumux {
namespace ShallowWater {

/*!
 * \brief Compute the outer cell state for fixed water depth boundary.
 *
 * \param waterDepthBoundary Discharge per meter at the boundary face [m^2/s]
 * \param waterDepthInside Water depth in the inner cell [m]
 * \param velocityXInside Velocity in x-direction in the inner cell [m/s]
 * \param velocityYInside Velocity in y-direction in the inner cell [m/s]
 * \param gravity Gravity constant [m/s^2]
 * \param nxy Normal vector of the boundary face
 *
 * \return cellStateOutside The outer cell state
 */
template<class Scalar, class GlobalPosition>
std::array<Scalar, 3> fixedWaterDepthBoundary(const Scalar waterDepthBoundary,
                                              const Scalar waterDepthInside,
                                              const Scalar velocityXInside,
                                              const Scalar velocityYInside,
                                              const Scalar gravity,
                                              const GlobalPosition& nxy)

{
    std::array<Scalar, 3> cellStateOutside;
    cellStateOutside[0] = waterDepthBoundary;

    using std::sqrt;
    const auto uboundIn = nxy[0] * velocityXInside  + nxy[1] * velocityYInside;
    const auto uboundQut =  uboundIn + 2.0 * sqrt(gravity * waterDepthInside) - 2.0 * sqrt(gravity * cellStateOutside[0]);

    cellStateOutside[1] = (nxy[0] * uboundQut); // we only use the normal part
    cellStateOutside[2] = (nxy[1] * uboundQut); // we only use the normal part

    return cellStateOutside;
}

/*!
 * \brief Compute the outer cell state for a fixed discharge boundary.
 *
 * \param dischargeBoundary Discharge per meter at the boundary face [m^2/s]
 * \param waterDepthInside Water depth in the inner cell [m]
 * \param velocityXInside Velocity in x-direction in the inner cell [m/s]
 * \param velocityYInside Velocity in y-direction in the inner cell [m/s]
 * \param gravity Gravity constant [m/s^2]
 * \param nxy Normal vector of the boundary face
 *
 * \return cellStateOutside The outer cell state
 */
template<class Scalar, class GlobalPosition>
std::array<Scalar, 3> fixedDischargeBoundary(const Scalar dischargeBoundary,
                                             const Scalar waterDepthInside,
                                             const Scalar velocityXInside,
                                             const Scalar velocityYInside,
                                             const Scalar gravity,
                                             const GlobalPosition& nxy)
{
    std::array<Scalar, 3> cellStateOutside;
    using std::abs;
    using std::sqrt;
    using std::max;

    // only impose if abs(q) > 0
    if (abs(dischargeBoundary) > 1.0e-9)
    {
        const auto uboundIn = nxy[0]*velocityXInside + nxy[1]*velocityYInside;
        const auto alphal = uboundIn + 2.0*sqrt(gravity * waterDepthInside);

        //initial guess for hstar solved with newton
        constexpr Scalar tol_hstar = 1.0E-12;
        constexpr Scalar ink_hstar = 1.0E-9;
        constexpr int maxstep_hstar = 30;

        Scalar hstar = 0.1;
        for (int i = 0; i < maxstep_hstar; ++i)
        {
            Scalar f_hstar = alphal - dischargeBoundary/hstar - 2 * sqrt(gravity * hstar);
            Scalar df_hstar = (f_hstar -(alphal - dischargeBoundary/(hstar + ink_hstar) - 2 * sqrt(gravity * (hstar+ink_hstar))))/ink_hstar;
            Scalar dx_hstar = -f_hstar/df_hstar;
            hstar = max(hstar - dx_hstar,0.001);

            if (dx_hstar*dx_hstar < tol_hstar)
                break;
        }

        const auto qinner = (nxy[0] * waterDepthInside * velocityYInside) - (nxy[1] * waterDepthInside * velocityXInside);
        cellStateOutside[0] = hstar;
        cellStateOutside[1] = (nxy[0] * dischargeBoundary - nxy[1] * qinner)/hstar;
        cellStateOutside[2] = (nxy[1] * dischargeBoundary + nxy[0] * qinner)/hstar;
    }

    return cellStateOutside;
}

/*!
 * \brief Compute the viscosity/diffusive flux at a rough wall boundary using no-slip formulation.
 *
 * \param alphaWall Roughness parameter: alphaWall=0.0 means full slip, alphaWall=1.0 means no slip, 0.0<alphaWall<1.0 means partial slip [-]
 * \param waterDepthInside Water depth in the inner cell [m]
 * \param velocityXInside Velocity in x-direction in the inner cell [m/s]
 * \param velocityYInside Velocity in y-direction in the inner cell [m/s]
 * \param turbulentViscosity Turbulent viscosity [m^2/s]
 * \param distance Distance from inside cell center to the boundary face [m]
 * \param nxy Normal vector of the boundary face
 *
 * \return cellStateOutside The outer cell state
 */
template<class Scalar, class GlobalPosition>
std::array<Scalar, 3> noslipWallBoundary(const Scalar alphaWall,
                                         const Scalar waterDepthInside,
                                         const Scalar velocityXInside,
                                         const Scalar velocityYInside,
                                         const Scalar turbulentViscosity,
                                         const Scalar distance,
                                         const GlobalPosition& nxy)
{
    std::array<Scalar, 3> roughWallFlux;
    using std::abs;

    // only impose if abs(alphaWall) > 0
    if (abs(alphaWall) > 1.0e-9)
    {
        // Initialization
        Scalar gradU   = 0.0;
        Scalar gradV   = 0.0;

        // Change of distance to wall for boundary conditions from Jasak (1996), p. 93-94:
        //Scalar dn[2] = {0.0};
        //Scalar fac = nxy[0]*dx+nxy[1]*dy;
        //dn[0] = fac*nxy[0];
        //dn[1] = fac*nxy[1];
        //distance = std::sqrt(dn[0]*dn[0] + dn[1]*dn[1]);

        // Compute the velocity gradients
        // Outside - inside cell: therefore the minus-sign
        // Only when cell contains sufficient water.
        // Use LET-limiter instead for differentiability?
        if ((waterDepthInside > 0.001))
        {
            gradU         = -alphaWall * velocityXInside/distance;
            gradV         = -alphaWall * velocityYInside/distance;
        }

        //At walls we assume the connection between the two cell centres to be
        //orthogonal to the boundary face, i.e. c_delta = 1.0
        Scalar c_delta = 1.0/(dx*nxy[0] + dy*nxy[1]);

        // Compute the viscosity/diffusive fluxes at the rough wall
        roughWallFlux[0] = 0.0;
        roughWallFlux[1] = turbulentViscosity * waterDepthInside * c_delta * gradU;
        roughWallFlux[2] = turbulentViscosity * waterDepthInside * c_delta * gradV;

    return roughWallFlux;
}

/*!
 * \brief Compute the viscosity/diffusive flux at a rough wall boundary using Nikuradse formulation.
 *
 * \param ksWall Nikuradse roughness height for the wall [m]
 * \param waterDepthInside Water depth in the inner cell [m]
 * \param velocityXInside Velocity in x-direction in the inner cell [m/s]
 * \param velocityYInside Velocity in y-direction in the inner cell [m/s]
 * \param distance Distance from inside cell center to the boundary face [m]
 * \param nxy Normal vector of the boundary face
 *
 * \return cellStateOutside The outer cell state
 */
template<class Scalar, class GlobalPosition>
std::array<Scalar, 3> nikuradseWallBoundary(const Scalar ksWall,
                                            const Scalar waterDepthInside,
                                            const Scalar velocityXInside,
                                            const Scalar velocityYInside,
                                            const Scalar distance,
                                            const GlobalPosition& nxy)
{
    std::array<Scalar, 3> roughWallFlux;
    using std::abs;
    using std::sqrt;
    using std::max;
    using std::log;

    // only impose if abs(ksWall) > 0
    if (abs(ksWall) > 1.0e-9)
    {
        Scalar y0w = ksWall/30.0;
        Scalar kappa2 = 0.41*0.41;
        // velocity magnitude
        auto velocityMagnitude  = sqrt(velocityXInside*velocityXInside + velocityYInside*velocityYInside);
        // and the wall shear stress components
        // should distance/y0w be limited to not become too small?
        auto tauWx = kappa2*velocityMagnitude*velocityXInside/ max(0.01,(log(distance/y0w + 1.0))*(log(distance/y0w + 1.0)));
        auto tauWy = kappa2*velocityMagnitude*velocityYInside/ max(0.01,(log(distance/y0w + 1.0))*(log(distance/y0w + 1.0)));

        // Compute the viscosity/diffusive fluxes at the rough wall
        roughWallFlux[0] = 0.0;
        roughWallFlux[1] = -waterDepthInside*tauWx; //*velocityXInside/max(0.01,abs(velocityXInside));
        roughWallFlux[2] = -waterDepthInside*tauWy; //*velocityYInside/max(0.01,abs(velocityYInside));

    return roughWallFlux;
}

} // end namespace ShallowWater
} // end namespace Dumux

#endif
