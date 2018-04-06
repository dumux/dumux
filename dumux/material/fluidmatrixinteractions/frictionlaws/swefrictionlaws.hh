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
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of different friction laws for the shallow water
 *        equations. This class bundles all friction laws into one file, in
 *        practical applications, the friction value(s) are stored at the
 *        nodes/cell centers. One can also use different friction laws in
 *        one application.
 *
 */
#ifndef DUMUX_SWE_FRICTIONLAWS_HH
#define DUMUX_SWE_FRICTIONLAWS_HH
#include <algorithm>
#include <cmath>

namespace Dumux
{
/*!
 * \ingroup Fluidmatrixinteractions
 *
 * \brief Implementation of different friction laws for the shallow water
 *        equations.
 *
 *
 */
template <class RF>
RF computeUstarH(const RF ks, const RF h, const RF grav, int law)
{

    /**************************************************************
    /
    / Friction as function of the water depth/gradient might work
    / better see
    /
    / Heniche et al. 2000 had a similar approach
    / "A two-dimensional finite element drying-wetting shallow
    / water model for rives and estuaries",
    / Advances in Water Resources 23(4):359-372
    /
    / We define a water depth minUpperH, if the water depth is
    / smaller, we start to limit the friciton law.
    / So the friciton term get's not extreme large for small water
    / depths
    /
    / ------------------------- minUpperh -----------
    /
    /
    /
    /------------------------rough_h ---------------
    /   /\  /\   roughness                  /grain\
    /-------------------------------bottom ------------------
    //////////////////////////////////////////////////
    /
    **************************************************************/

    RF ustar_h, rough_h;
    using std::pow;
    using std::max;
    using std::min;
    using std::log;

    const int manning = 1;
    const int chezy = 2;
    const int nikuradse = 3;

    const RF letL = 0.0, letT = 2.0, letE = 1.0; //let this parameters as they are
    RF mobility = 1.0;
    RF krw = 1.0;
    RF sw = 0.0;

    ustar_h = 0.0;    //compute friction


    if (law == manning){
        rough_h = std::pow(25.68/(1.0/ks),6.0);
    }else{
        rough_h = ks;
    }

    auto minUpperH = rough_h* 2.0; //we start to limit the friction for small water depths

    //Mobility after LET-model
    sw = min(h * (1.0/minUpperH),1.0);
    sw = max(0.0,sw);
    mobility = (krw * pow(sw,letL))/(pow(sw,letL) + letE * pow(1.0-sw,letT));
    rough_h = rough_h * (1.0 - mobility);

    //compute friction
    if (law == manning){
        //Manning friction
        if (ks > 0){
            auto cfric = pow((h + rough_h),1.0/6.0) * 1.0/(ks);
            ustar_h = grav / pow(cfric,2.0);
        }
    }else if(law == chezy){
        //Chezy friction law
        ustar_h = grav / pow(ks,2.0);

    }else if(law == nikuradse){
        //Nikuradse law
        ustar_h = pow(0.41,2.0)/pow(log((12*(h+ rough_h))/ks),2.0);
    }else{
        //now law is used and ustar_h is zero
        ustar_h = 0.0;
    }

    return ustar_h;
}

template <class RF>
RF computeKsH(const RF ks, int law)
{

    /**************************************************************
    / Compute an estimated height of the friction elements
    / (grains, etc.)
    /
    / This is a numerical parameter to limit the flux for depths
    / smaller as the estimated height
    **************************************************************/

    RF ks_h;
    using std::pow;
    const int manning = 1;
    const int chezy = 2;
    const int nikuradse = 3;

    if (law == manning){
        ks_h = std::pow(25.68/(1.0/ks),6.0);
    }else{
        ks_h = ks;
    }

    return ks_h;
}



} // end namespace Dumux

#endif // SWE_FRICTIONLAWS_HH
