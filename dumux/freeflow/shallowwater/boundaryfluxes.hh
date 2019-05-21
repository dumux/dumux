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
 * \brief Compute boundary conditions (cell state) via Riemann invariants
 *
 */
#ifndef DUMUX_SHALLOWWATER_BOUNDARYFLUXES_HH
#define DUMUX_SHALLOWWATER_BOUNDARYFLUXES_HH

#include <array>
#include <cmath>

namespace Dumux {
namespace ShallowWater {

/*!
 * \brief compute the cell state for fixed water depth boundary.
 */
template<class Scalar, class GlobalPosition>
std::array<Scalar, 3> fixedWaterDepthBoundary(const Scalar waterDepthBoundary,
                                            const Scalar waterDepthLeft,
                                            const Scalar velocityXLeft,
                                            const Scalar velocityYLeft,
                                            const GlobalPosition& nxy)

{
    std::array<Scalar, 3> cellStateRight;
    cellStateRight[0] = waterDepthBoundary;

    using std::sqrt;
    const auto uboundIn = nxy[0] * velocityXLeft  + nxy[1] * velocityYLeft;
    const auto uboundQut =  uboundIn + 2.0 * sqrt(9.81 * waterDepthLeft) - 2.0 * sqrt(9.81 * cellStateRight[0]);

    cellStateRight[1] = (nxy[0] * uboundQut); // we only use the normal part
    cellStateRight[2] = (nxy[1] * uboundQut); // we only use the normal part

    return cellStateRight;
}

/*!
 * \brief compute the cell state for a fixed discharge boundary.
 */
template<class Scalar, class GlobalPosition>
std::array<Scalar, 3> fixedDischargeBoundary(const Scalar qlocal,
                                           const Scalar waterDepthLeft,
                                           const Scalar velocityXLeft,
                                           const Scalar velocityYLeft,
                                           const GlobalPosition& nxy)
{
    std::array<Scalar, 3> cellStateRight;
    using std::abs;
    using std::sqrt;
    using std::max;

    // only impose if abs(q) > 0
    if (abs(qlocal) > 1.0e-9)
    {
        const auto uboundIn = nxy[0]*velocityXLeft + nxy[1]*velocityYLeft;
        const auto alphal = uboundIn + 2.0*sqrt(9.81 * waterDepthLeft);

        //initial guess for hstar solved with newton
        constexpr Scalar tol_hstar = 1.0E-12;
        constexpr Scalar ink_hstar = 1.0E-9;
        constexpr int maxstep_hstar = 30;

        Scalar hstar = 0.1;
        for (int i = 0; i < maxstep_hstar; ++i)
        {
            Scalar f_hstar = alphal - qlocal/hstar - 2 * sqrt(9.81 * hstar);
            Scalar df_hstar = (f_hstar -(alphal - qlocal/(hstar + ink_hstar) - 2 * sqrt(9.81 * (hstar+ink_hstar))))/ink_hstar;
            Scalar dx_hstar = -f_hstar/df_hstar;
            hstar = max(hstar - dx_hstar,0.001);

            if (dx_hstar*dx_hstar < tol_hstar)
                break;
        }

        const auto qinner = (nxy[0] * waterDepthLeft * velocityYLeft) - (nxy[1] * waterDepthLeft * velocityXLeft);
        cellStateRight[0] = hstar;
        cellStateRight[1] = (nxy[0] * qlocal - nxy[1] * qinner)/hstar;
        cellStateRight[2] = (nxy[1] * qlocal + nxy[0] * qinner)/hstar;
    }

    return cellStateRight;
}

} // end namespace ShallowWater
} // end namespace Dumux

#endif
