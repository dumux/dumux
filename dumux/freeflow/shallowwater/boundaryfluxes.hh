// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ShallowWaterModels
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
#include <algorithm>

namespace Dumux::ShallowWater {

/*!
 * \brief Compute the outer cell state for a slip wall.
 *
 * The normal velocity is mirrored and the tangential velocity kept, so the Riemann
 * problem across the face carries no mass flux and the hydrostatic pressure force.
 *
 * \param waterDepthInside Water depth in the inner cell [m]
 * \param velocityXInside Velocity in x-direction in the inner cell [m/s]
 * \param velocityYInside Velocity in y-direction in the inner cell [m/s]
 * \param nxy Normal vector of the boundary face
 *
 * \return cellStateOutside The outer cell state
 */
template<class Scalar, class GlobalPosition>
std::array<Scalar, 3> wallBoundary(const Scalar waterDepthInside,
                                   const Scalar velocityXInside,
                                   const Scalar velocityYInside,
                                   const GlobalPosition& nxy)
{
    const auto uNormalIn = nxy[0]*velocityXInside + nxy[1]*velocityYInside;
    return {{waterDepthInside,
             velocityXInside - 2.0*uNormalIn*nxy[0],
             velocityYInside - 2.0*uNormalIn*nxy[1]}};
}

/*!
 * \brief Compute the outer cell state for an inflow with prescribed depth and speed.
 *
 * Both characteristics enter the domain when the inflow is supercritical, so depth
 * and velocity have to be imposed together. For subcritical inflow one of them is
 * determined by the interior, see fixedDischargeBoundary and fixedWaterDepthBoundary.
 *
 * \param waterDepthBoundary Water depth at the boundary face [m]
 * \param speedBoundary Speed of the inflowing water, positive into the domain [m/s]
 * \param nxy Normal vector of the boundary face
 *
 * \return cellStateOutside The outer cell state
 */
template<class Scalar, class GlobalPosition>
std::array<Scalar, 3> inflowBoundary(const Scalar waterDepthBoundary,
                                     const Scalar speedBoundary,
                                     const GlobalPosition& nxy)
{
    return {{waterDepthBoundary, -speedBoundary*nxy[0], -speedBoundary*nxy[1]}};
}

/*!
 * \brief Compute the outer cell state for a free overfall (open outflow).
 *
 * The boundary characteristics travel with the normal velocity plus or minus the wave
 * celerity. Supercritical outflow carries no characteristic back into the domain, so the
 * outer state is the inner one. Subcritical outflow does carry one; reusing the inner state
 * then imposes a zero gradient, which under-drains and piles water up in the boundary
 * cell, whereas the brink of a free overfall sits at critical depth for the arriving
 * normal discharge. The tangential velocity is kept. Water cannot flow up the brink,
 * so a non-positive normal velocity gives the wall state.
 *
 * \param waterDepthInside Water depth in the inner cell [m]
 * \param velocityXInside Velocity in x-direction in the inner cell [m/s]
 * \param velocityYInside Velocity in y-direction in the inner cell [m/s]
 * \param gravity Gravity constant [m/s^2]
 * \param nxy Normal vector of the boundary face
 *
 * \return cellStateOutside The outer cell state
 */
template<class Scalar, class GlobalPosition>
std::array<Scalar, 3> criticalDepthOutflowBoundary(const Scalar waterDepthInside,
                                                   const Scalar velocityXInside,
                                                   const Scalar velocityYInside,
                                                   const Scalar gravity,
                                                   const GlobalPosition& nxy)
{
    using std::sqrt; using std::cbrt; using std::max;

    const auto uNormalIn = nxy[0]*velocityXInside + nxy[1]*velocityYInside;
    if (uNormalIn <= 0.0)
        return wallBoundary(waterDepthInside, velocityXInside, velocityYInside, nxy);

    const auto froude = uNormalIn/sqrt(gravity*max(waterDepthInside, Scalar(1e-12)));
    if (froude >= 1.0)
        return {{waterDepthInside, velocityXInside, velocityYInside}};

    const auto normalDischarge = waterDepthInside*uNormalIn;
    const auto criticalDepth = cbrt(normalDischarge*normalDischarge/gravity);
    const auto uNormalOut = sqrt(gravity*criticalDepth);
    const auto uTangentialIn = -nxy[1]*velocityXInside + nxy[0]*velocityYInside;
    return {{criticalDepth,
             nxy[0]*uNormalOut - nxy[1]*uTangentialIn,
             nxy[1]*uNormalOut + nxy[0]*uTangentialIn}};
}

/*!
 * \brief Compute the outer cell state for fixed water depth boundary.
 *
 * \param waterDepthBoundary Water depth at the boundary face [m^2/s]
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
    const auto uNormalIn = nxy[0] * velocityXInside + nxy[1] * velocityYInside;
    const auto uTangentialIn = -nxy[1] * velocityXInside + nxy[0] * velocityYInside;
    const auto uNormalOut = uNormalIn + 2.0 * sqrt(gravity * waterDepthInside) - 2.0 * sqrt(gravity * cellStateOutside[0]);

    cellStateOutside[1] = (nxy[0] * uNormalOut - nxy[1] * uTangentialIn);
    cellStateOutside[2] = (nxy[1] * uNormalOut + nxy[0] * uTangentialIn);

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
    using std::abs;
    using std::sqrt;
    using std::max;

    // a vanishing discharge is a wall
    auto cellStateOutside = wallBoundary(waterDepthInside, velocityXInside, velocityYInside, nxy);

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
 * \param alphaWall Roughness parameter: alphaWall=0.0 means full slip, alphaWall=1.0 means no slip, \f$0.0<\text{alphaWall}<1.0\f$ means partial slip [-]
 * \param turbulentViscosity Turbulent viscosity [m^2/s]
 * \param state Primary variables (water depth, velocities)
 * \param cellCenterToBoundaryFaceCenter Cell-center to boundary distance
 * \param unitNormal Normal vector of the boundary face
 */
template<class PrimaryVariables, class Scalar, class GlobalPosition>
std::array<Scalar, 3> noslipWallBoundary(const Scalar alphaWall,
                                         const Scalar turbulentViscosity,
                                         const PrimaryVariables& state,
                                         const GlobalPosition& cellCenterToBoundaryFaceCenter,
                                         const GlobalPosition& unitNormal)
{
    // only impose if abs(alphaWall) > 0
    using std::abs;
    if (abs(alphaWall) <= 1.0e-9)
        return {};

    const auto waterDepth = state[0];
    // regularization: Set gradients to zero for drying cell
    // Use LET-limiter instead for differentiability?
    if (waterDepth <= 0.001)
        return {};

    const auto xVelocity  = state[1];
    const auto yVelocity  = state[2];
    const auto distance = cellCenterToBoundaryFaceCenter.two_norm();

    // Compute the velocity gradients
    // Outside - inside cell: therefore the minus-sign
    // Only when cell contains sufficient water.
    const auto gradU = -alphaWall * xVelocity/distance;
    const auto gradV = -alphaWall * yVelocity/distance;

    // Factor that takes the direction of the unit vector into account
    const auto direction = (unitNormal*cellCenterToBoundaryFaceCenter)/distance;

    // Compute the viscosity/diffusive fluxes at the rough wall
    return {
        0.0,
        -turbulentViscosity*waterDepth*gradU*direction,
        -turbulentViscosity*waterDepth*gradV*direction
    };
}

/*!
 * \brief Compute the viscosity/diffusive flux at a rough wall boundary using Nikuradse formulation.
 *
 * \param ksWall Nikuradse roughness height for the wall [m]
 * \param state the primary variable state (water depth, velocities)
 * \param cellCenterToBoundaryFaceCenter Cell-center to boundary distance
 * \param unitNormal Normal vector of the boundary face
 */
template<class PrimaryVariables, class Scalar, class GlobalPosition>
std::array<Scalar, 3> nikuradseWallBoundary(const Scalar ksWall,
                                            const PrimaryVariables& state,
                                            const GlobalPosition& cellCenterToBoundaryFaceCenter,
                                            const GlobalPosition& unitNormal)
{
    // only impose if abs(ksWall) > 0
    using std::abs;
    if (abs(ksWall) <= 1.0e-9)
        return {};

    using std::hypot;
    const Scalar xVelocity = state[1];
    const Scalar yVelocity = state[2];
    const Scalar velocityMagnitude = hypot(xVelocity, yVelocity);
    const Scalar distance = cellCenterToBoundaryFaceCenter.two_norm();
    const Scalar y0w = ksWall/30.0;
    constexpr Scalar kappa2 = 0.41*0.41;

    // should distance/y0w be limited to not become too small?
    using std::log; using std::max;
    const auto logYPlus = log(distance/y0w+1.0);
    const auto fac = kappa2*velocityMagnitude / max(1.0e-3,logYPlus*logYPlus);

    // Factor that takes the direction of the unit vector into account
    const auto direction = (unitNormal*cellCenterToBoundaryFaceCenter)/distance;

    // wall shear stress vector
    const auto tauWx = direction*fac*xVelocity;
    const auto tauWy = direction*fac*yVelocity;

    // Compute the viscosity/diffusive fluxes at the rough wall
    const auto waterDepth = state[0];
    return {0.0, waterDepth*tauWx, waterDepth*tauWy};
}

} // end namespace Dumux::ShallowWater

#endif
