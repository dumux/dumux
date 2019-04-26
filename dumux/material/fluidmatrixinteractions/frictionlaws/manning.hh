// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \brief Implementation of the friction law after Manning.
 *
 * The LET mobility model is used to limit the friction for small water depths.
 */
#ifndef DUMUX_FRICTIONLAW_MANNING_HH
#define DUMUX_FRICTIONLAW_MANNING_HH

#include <algorithm>
#include <cmath>

namespace Dumux {
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the friction law after Manning.
 *
 * The LET mobility model is used to limit the friction for small water depths.
 */

template <typename RF>
class FrictionLawManning
{
public:
    /*!
     * \brief  The friction law after Manning.
     * \param volVars volume variables
     */

    template<class VolumeVariables>
    RF computeUstarH(const VolumeVariables& volVars, RF ks)
    {
        using std::pow;
        using std::min;
        using std::max
        RF ustar_h = 0.0;
        RF rough_h = pow(25.68/(1.0/ks),6.0);

        // Limit the friction for small water depths with the LET model.
        // ???Proberply timo will complain about the names???
        const RF letL = 0.0;
        const RF letT = 2.0;
        const RF letE = 1.0;
        RF mobility = 1.0;
        RF krw = 1.0;
        RF sw = 0.0;
        auto minUpperH = rough_h * 2.0;
        sw = min(volVars::getH * (1.0/minUpperH),1.0);
        sw = max(0.0,sw);
        mobility = (krw * pow(sw,letL))/(pow(sw,letL) + letE * pow(1.0-sw,letT));
        rough_h = rough_h * (1.0 - mobility);

        auto cfric = pow((volVars::getH + rough_h),1.0/6.0) * 1.0/(ks);
        ustar_h = volVars::getGravity / pow(cfric,2.0);
        return ustar_h;
    }
};

} // end namespace Dumux

#endif // DUMUX_FRICTIONLAW_MANNING_HH
