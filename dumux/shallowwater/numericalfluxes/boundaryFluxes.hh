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
 * \ingroup SweModel
 * \brief Compute the boundary fluxes based on the Riemann invariants
 *
 *        bdType can be:
 *          0 = noFlow (do nothing)
 *          1 = given discharge
 *          2/3 = given h (if 3, h comes from a stage discharge curve)
 *
 *
 */
#ifndef DUMUX_SHALLOWWATER_NUMERICALFLUXES_BOUNDARYFLUXES_HH
#define DUMUX_SHALLOWWATER_NUMERICALFLUXES_BOUNDARYFLUXES_HH


#include <dumux/common/math.hh>

namespace Dumux
{

inline void boundaryFluxes(const auto& nxy,
                               const auto& faceVolume,
                               const auto& minH,
                               const auto& cellStatesLeft,
                               auto& cellStatesRight,
                               const auto& bdType,
                               const auto& bdValue)
  {
    if (bdType==0)//no flow
    {
        //change nothing since it is the default case
    }
    //*********         set given discharge                     **********//
    if (bdType==1)//defined q
    {
      //olny impose if abs(q) > 0
      if (std::abs(bdValue) > 1.0E-9){
        auto qlocal =  (bdValue) /faceVolume;
        auto uboundIn = nxy[0] * cellStatesLeft[1] + nxy[1] * cellStatesLeft[2];
        auto alphal = uboundIn + 2.0 * std::sqrt(9.81 * cellStatesLeft[0]);

        //initial guess for hstar solved with newton
        double hstar = 0.1;
        double tol_hstar = 1.0E-12;
        double f_hstar;
        double df_hstar;
        double ink_hstar = 1.0E-9;
        double dx_hstar;
        int maxstep_hstar = 30;

        for(int i = 0; i < maxstep_hstar; ++i){
            f_hstar = alphal - qlocal/hstar - 2*std::sqrt(9.81 * hstar);
            df_hstar = (f_hstar -(alphal - qlocal/(hstar + ink_hstar) - 2*std::sqrt(9.81 * (hstar+ink_hstar))))/ink_hstar;
            dx_hstar = -f_hstar/df_hstar;
            hstar = std::max(hstar - dx_hstar,0.001);

            if (std::pow(dx_hstar,2.0) < tol_hstar){
                    break;
            }
         }
         auto qinner = (nxy[0] * cellStatesLeft[0] * cellStatesLeft[2]) - (nxy[1] * cellStatesLeft[0] * cellStatesLeft[1]);
         cellStatesRight[0] = hstar;
         cellStatesRight[1] = (nxy[0] * qlocal - nxy[1] * qinner)/hstar;
         cellStatesRight[2] = (nxy[1] * qlocal + nxy[0] * qinner)/hstar;
      }
    }
    //*********  set fixed h **********//
    if (((bdType==2)||(bdType==3))&&(bdValue - cellStatesRight[3] > minH))
    {
      cellStatesRight[0] = bdValue - cellStatesRight[3];
      auto uboundIn = nxy[0] * cellStatesLeft[1]  + nxy[1] * cellStatesLeft[2] ;
      auto uboundQut =  uboundIn + 2.0 * std::sqrt(9.81 * cellStatesLeft[0]) - 2.0 * std::sqrt(9.81 * cellStatesRight[0]);
      cellStatesRight[1] = (nxy[0] * uboundQut); // - nxy[1] * (-nxy[1] * cellStatesLeft[1] + nxy[0] * cellStatesLeft[1]));
      cellStatesRight[2] = (nxy[1] * uboundQut); // + nxy[0] * (-nxy[1] * cellStatesLeft[1] + nxy[0] * cellStatesLeft[1]));
    }
  }
} // end namespace Dumux

#endif
