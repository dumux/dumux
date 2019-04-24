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
 * \ingroup ShallowWater
 * \brief Compute the boundary fluxes based on the Riemann invariants
 *
 */
#ifndef DUMUX_SHALLOWWATER_BOUNDARYFLUXES_HH
#define DUMUX_SHALLOWWATER_BOUNDARYFLUXES_HH

#include <array>
#include <algorithm>
#include <dumux/common/math.hh>

namespace Dumux {
namespace ShallowWater {

template<class Scalar, class GlobalPosition>
std::array<Scalar,3> fixedWaterDepthBoundary(Scalar waterDepthBoundary,
                                             Scalar waterDepthLeft,
                                             Scalar waterDepthRight,
                                             Scalar velocityXLeft,
                                             Scalar velocityXRight,
                                             Scalar velocityYLeft,
                                             Scalar velocityYRight,
                                             Scalar gravity,
                                             GlobalPosition nxy)

{
    std::array<Scalar,3> cellStateRight;

    cellStateRight[0] = waterDepthBoundary;

    auto uboundIn = nxy[0] * velocityXLeft  + nxy[1] * velocityYLeft ;
    auto uboundQut =  uboundIn + 2.0 * sqrt(9.81 * waterDepthLeft) - 2.0 * sqrt(9.81 * cellStateRight[0]);
    cellStateRight[1] = (nxy[0] * uboundQut); // we only use the normal part
    cellStateRight[2] = (nxy[1] * uboundQut); // we only use the normal part

    return cellStateRight;
}

template<class Scalar, class GlobalPosition>
std::array<Scalar,3> fixedDischargeBoundary(Scalar dischargeBoundary,
                                            Scalar waterDepthLeft,
                                            Scalar waterDepthRight,
                                            Scalar velocityXLeft,
                                            Scalar velocityXRight,
                                            Scalar velocityYLeft,
                                            Scalar velocityYRight,
                                            Scalar gravity,
                                            GlobalPosition nxy,
                                            Scalar faceVolume)
{
    std::array<Scalar,3> cellStateRight;
    using std::pow;
    using std::abs;
    using std::sqrt;

    //olny impose if abs(q) > 0
    if (abs(dischargeBoundary) > 1.0E-9){
        auto qlocal =  (dischargeBoundary) /faceVolume;
        auto uboundIn = nxy[0] * velocityXLeft + nxy[1] * velocityYLeft;
        auto alphal = uboundIn + 2.0 * sqrt(9.81 * waterDepthLeft);

        //initial guess for hstar solved with newton
        Scalar hstar = 0.1;
        Scalar tol_hstar = 1.0E-12;
        Scalar ink_hstar = 1.0E-9;
        int maxstep_hstar = 30;

        for(int i = 0; i < maxstep_hstar; ++i){
            Scalar f_hstar = alphal - qlocal/hstar - 2 * sqrt(9.81 * hstar);
            Scalar df_hstar = (f_hstar -(alphal - qlocal/(hstar + ink_hstar) - 2 * sqrt(9.81 * (hstar+ink_hstar))))/ink_hstar;
            Scalar dx_hstar = -f_hstar/df_hstar;
            hstar = max(hstar - dx_hstar,0.001);

            if (pow(dx_hstar,2.0) < tol_hstar){
                break;
            }
        }
        auto qinner = (nxy[0] * waterDepthLeft * velocityYLeft) - (nxy[1] * waterDepthLeft * velocityXLeft);
        cellStateRight[0] = hstar;
        cellStateRight[1] = (nxy[0] * qlocal - nxy[1] * qinner)/hstar;
        cellStateRight[2] = (nxy[1] * qlocal + nxy[0] * qinner)/hstar;
    }

    return cellStateRight;
}

} // end namespace ShallowWater
} // end namespace Dumux

#endif
