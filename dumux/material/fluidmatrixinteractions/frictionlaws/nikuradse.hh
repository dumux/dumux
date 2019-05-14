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
 * \brief Implementation of the friction law after Nikuradse.
 *
 * The LET mobility model is used to limit the friction for small water depths.
 */
#ifndef DUMUX_FRICTIONLAW_NIKURADSE_HH
#define DUMUX_FRICTIONLAW_NIKURADSE_HH

#include <algorithm>
#include <cmath>

namespace Dumux {
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the friction law after Nikuradse.
 *
 * The LET mobility model is used to limit the friction for small water depths.
 */

template <typename Scalar>
class FrictionLawNikuradse
{
public:
    /*!
     * \brief Compute the friction ustar_h.
     *
     * \param h water depth.
     * \param ks the Strickler friction value.
     * \return ustar_h friction used for the source term in shallow water models.
     */

    Scalar computeUstarH(const Scalar h,const Scalar ks)
    {
        using std::pow;
        using std::log;

        Scalar ustar_h = 0.0;
        Scalar rough_h = ks;

        rough_h = limitRoughH(rough_h, h);
        ustar_h = pow(0.41,2.0)/pow(log((12*(h + rough_h))/ks),2.0);
        return ustar_h;
    }

    /*!
     * \brief Limit the friction for small water depth.
     *
     * We define a water depth minUpperH. If the water depth is
     * smaller, we start to limit the friciton.
     * So the friciton term get's not extreme large for small water
     * depths.
     *
     * ------------------------- minUpperh -----------
     *
     *
     *
     * ------------------------rough_h ---------------
     *    /\  /\   roughness                  /grain\
     * -------------------------------bottom ------------------
     * /////////////////////////////////////////////////
     *
     * For the limiting the LET model is used, which is usually applied in the
     * porous media flow to limit the permeability due to the saturation. It employs
     * the three empirical paramaters L, E and T, which describe the limiting curve.
     *
     * \param rough_h roughness height of the representive structure (e.g. largest grain size).
     * \param h water depth.
     */
    Scalar limitRoughH(Scalar rough_h, const Scalar h)
    {
        using std::pow;
        using std::min;
        using std::max;

        const Scalar letL = 0.0; //!< empirical parameter of the LET model
        const Scalar letT = 2.0; //!< empirical parameter of the LET model
        const Scalar letE = 1.0; //!< empirical parameter of the LET model
        Scalar mobility_max = 1.0; //!< maximal mobility

        auto minUpperH = rough_h * 2.0;
        auto sw = min(h * (1.0/minUpperH),1.0);
        sw = max(0.0,sw);
        auto mobility = (mobility_max * pow(sw,letL))/(pow(sw,letL) + letE * pow(1.0-sw,letT));
        return rough_h * (1.0 - mobility);
    }
};

} // end namespace Dumux

#endif // DUMUX_FRICTIONLAW_NIKURADSE_HH
