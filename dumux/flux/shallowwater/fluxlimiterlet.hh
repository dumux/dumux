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
 * \ingroup ShallowWaterFlux
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
 * \ingroup ShallowWaterFlux
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
 * \param upperH Where to start the limit the function (mobility < 1)
 * \param lowerH Where the limit should reach zero (mobility < 0)
 */
template<class Scalar>
static Scalar fluxLimiterLET(const Scalar valueLeft,
                             const Scalar valueRight,
                             const Scalar upperH,
                             const Scalar lowerH)
{
    using std::pow;
    using std::min;
    using std::max;

    const auto h = (valueLeft + valueRight)*0.5;

    Scalar mobility = 1.0;
    if (h < upperH)
    {
        const auto sw = max(min(h*(1.0/upperH) - lowerH, 1.0), 0.0);

        // LET-model for mobility
        // constexpr Scalar krw = 1.0;
        // constexpr Scalar letL = 2.0;
        // constexpr Scalar letT = 2.0;
        // constexpr Scalar letE = 1.0;
        // mobility = (krw * pow(sw, letL))/(pow(sw, letL) + letE * pow(1.0 - sw, letT));

        mobility = (sw*sw)/(sw*sw + (1-sw)*(1-sw));
    }

    return mobility;
}

} // end namespace ShallowWater
} // end namespace Dumux

#endif
