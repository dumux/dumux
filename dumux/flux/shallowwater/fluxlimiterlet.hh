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
 * \brief Function to limit the fluxes
 *
 */

#ifndef DUMUX_FLUX_SHALLOW_WATER_FLUX_LIMITER_LET_HH
#define DUMUX_FLUX_SHALLOW_WATER_FLUX_LIMITER_LET_HH

#include <algorithm>
#include <cmath>

namespace Dumux {
namespace ShallowWater {
/*!
 * \ingroup ShallowWater
 * \brief Flux limiter function to scale fluxes for small water depths.
 *
 * This function acts like a kind of mobility, it limits the water flux
 * for small water depths. The mobility depends on the left and right
 * side state of a variable. The LET-Type function is described at
 * https://en.wikipedia.org/wiki/Relative_permeability
 * The LET-Parameters are fixed. The current parameters are set to a
 * values to get a nice curve. They have no physical meaning.
 *
 * \tparam Scalar the scalar type for scalar physical quantities
 * \param valueLeft The value on the left side
 * \param valueRight The value on the right side
 * \param upperLimit Where to start the limit the function (mobility < 1)
 * \param lowerLimit Where the limit should reach zero (mobility < 0)
 */
template<class Scalar>
static Scalar fluxLimiterLET(const Scalar& valueLeft,
                             const Scalar& valueRight,
                             const Scalar& upperH,
                             const Scalar& lowerH)
{

        Scalar krw = 1.0;
        Scalar sw = 0.0;
        Scalar letL = 2.0;
        Scalar letT = 2.0;
        Scalar letE = 1.0;
        Scalar h = 0.0;
        Scalar mobility = 1.0;
        using std::pow;
        using std::min;
        using std::max;

        h = (valueLeft+valueRight)*0.5;

        if (h < upperH)
        {
            sw = min(h * (1.0/upperH) - (lowerH),1.0);
            sw = max(sw,0.0);
            sw = min(sw,1.0);

            //LET-model for mobility
            mobility = (krw * pow(sw, letL))/(pow(sw, letL) + letE * pow(1.0 - sw, letT));
        }

        return mobility;
}
} // end namespace ShallowWater
} // end namespace Dumux

#endif
