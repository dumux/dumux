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
 *
 * \brief Specification of threshold capillary pressures for the PNM.
 */
#ifndef DUMUX_PNM_THRESHOLD_CAPILLARY_PRESSURES_HH
#define DUMUX_PNM_THRESHOLD_CAPILLARY_PRESSURES_HH

#include <cmath>

namespace Dumux
{

class ThresholdCapillaryPressures
{
public:

    //! The snap-off capillary pressure of a pore throat
    template<class Scalar>
    static constexpr Scalar pcSnapoff(const Scalar surfaceTension,
                                      const Scalar contactAngle,
                                      const Scalar inscribedRadius) noexcept
    {
        using std::sin;
        using std::cos;
        const Scalar theta = contactAngle;
        const Scalar cosTheta = std::cos(theta);
        const Scalar sinTheta = std::sin(theta);
        return surfaceTension / inscribedRadius * (cosTheta - sinTheta);
    }

    /*! \brief The entry capillary pressure of a pore throat.
     *
     * For details, see Eq. 11 in Rabbani et al., 2016
     * or Eq A-7 in Oren et al., 1998
     */
    template<class Scalar>
    static constexpr Scalar pcEntry(const Scalar surfaceTension,
                                    const Scalar contactAngle,
                                    const Scalar inscribedRadius,
                                    const Scalar shapeFactor) noexcept
    {
        using std::sin;
        using std::cos;
        using std::sqrt;
        const Scalar theta = contactAngle;
        const Scalar cosTheta = cos(theta);
        const Scalar sinTheta = sin(theta);

        const Scalar D = M_PI - 3*theta + 3*sinTheta*cosTheta - (cosTheta*cosTheta) / (4*shapeFactor);
        const Scalar F = (1 + sqrt(1 + 4*shapeFactor*D / (cosTheta*cosTheta))) / (1 + 2*sqrt(M_PI*shapeFactor));
        return surfaceTension / inscribedRadius * cosTheta * (1 + 2*sqrt(M_PI*shapeFactor)) * F;
    }

};

} // namespace Dumux

#endif
