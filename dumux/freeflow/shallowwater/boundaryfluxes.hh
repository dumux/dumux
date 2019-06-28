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

} // end namespace ShallowWater
} // end namespace Dumux

#endif
