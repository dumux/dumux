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
 * \ingroup CCTpfaDiscretization
 * \brief Simple model to compute a mobility of water
 *        for small water depths the mobility gets zero.
 *        No limitation for the momentum terms!
 *        Parameters letL, letE, letT describe the curve, see
 *        see https://en.wikipedia.org/wiki/Relative_permeability.
 *
 *        For water fluxes we have dh as the actual value, for water
 *        we use (dl+dr)/2.0,
 *        dhmin = ks * 0.01;
 *        dhmax = ks * 1.0;
 *
 *
 */
#ifndef DUMUX_SHALLOWWATER_NUMERICALFLUXES_LETMODEL_HH
#define DUMUX_SHALLOWWATER_NUMERICALFLUXES_LETMODEL_HH


#include <dumux/common/math.hh>

namespace Dumux
{


inline void letmobility(const auto dl,const auto dr, const auto ks, auto& mobility)
  {

        //Add a mobility term
        //LET-Type see https://en.wikipedia.org/wiki/Relative_permeability
        //This is hard coded so far
        double krw = 1.0;
        double sw = 0.0;
        double minUpperH,minLowerH;
        double letL = 2.0;
        double letT = 2.0;
        double letE = 1.0;
        mobility[0] = 1.0;
        mobility[1] = 1.0;
        mobility[2] = 1.0;

        minUpperH = ks*1.0; //We use ks*1 and compute the maximum in ks_av
        minLowerH = ks*0.01; //0.001 is a good lower limit

        double h;

        h = (dl+dr)*0.5;
        //h = std::min(dl,dr);
        //h = std::max(dl,dr);

        if (h < minUpperH){
          sw = std::min(h * (1.0/minUpperH) - (minLowerH),1.0);  //best factor is not totally clear, about 0.01?!
          sw = std::max(sw,0.0);
          sw = std::min(sw,1.0);

          //LET-model for mobility
          mobility[0] = (krw * std::pow(sw,letL))/(std::pow(sw,letL) + letE* std::pow(1.0-sw,letT));
        }

    }
} // end namespace Dumux

#endif
