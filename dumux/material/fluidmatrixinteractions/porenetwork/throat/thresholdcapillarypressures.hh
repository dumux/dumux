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
 * \ingroup PoreNetworkModels
 * \brief Specification of threshold capillary pressures for the PNM.
 */
#ifndef DUMUX_PNM_THRESHOLD_CAPILLARY_PRESSURES_HH
#define DUMUX_PNM_THRESHOLD_CAPILLARY_PRESSURES_HH

#include <cmath>

#include <dune/common/exceptions.hh>
#include <dumux/porenetwork/common/throatproperties.hh>

namespace Dumux::PoreNetwork::ThresholdCapillaryPressures {

/*! \brief The snap-off capillary pressure of a pore throat with regular cross section shape
 * (with same corner angles and side length)
 * For details, see Eq. 4.8 in Multiphase Flow in
 * Blunt, M. J. (2017). Multiphase flow in permeable media: A pore-scale perspective. Cambridge university press.
 * https://doi.org/10.1017/9781316145098
 */
template <class Scalar>
Scalar pcSnapoffRegularShape(const Scalar surfaceTension,
                             const Scalar contactAngle,
                             const Scalar inscribedRadius,
                             const Throat::Shape shape) noexcept
{
    using std::cos;
    using std::sin;

    const Scalar cornerHalfAngle = Throat::cornerHalfAngles<Scalar>(shape)[0];

    // check if snap-off is possible with such contact angle and corner half-angle
    // if the criteriun is not fulfilled, return the lowest negative value to make sure that snap-off will not happen
    if (cornerHalfAngle + contactAngle >= 0.5*M_PI)
        return std::numeric_limits<Scalar>::lowest();

    const Scalar cosContactAngle = cos(contactAngle);
    const Scalar sinContactAngle = sin(contactAngle);

    // Get tangent of corner half-angle
    const Scalar tanCornerHalfAngle = [&]{ switch (shape){ // optimized tan(corner half-angle) for certain shapes
        case Throat::Shape::equilateralTriangle: return 0.577;
        case Throat::Shape::square: return 1.0;
        default: using std::tan; return tan(cornerHalfAngle);
    }}();

    return surfaceTension / inscribedRadius * (cosContactAngle - sinContactAngle * tanCornerHalfAngle);
}


/*! \brief The snap-off capillary pressure of a pore throat
 * It checks if the cross section shape of the throat is a regular or irregular shape and call
 * the proper pc snap-off accordingly
 */
template <class Scalar>
Scalar pcSnapoff(const Scalar surfaceTension,
                 const Scalar contactAngle,
                 const Scalar inscribedRadius,
                 const Throat::Shape shape)
{
    if (Throat::isRegularShape(shape)) // call the pc snap-off calculated for regular shapes(same corner angles and side length)
        return pcSnapoffRegularShape(surfaceTension, contactAngle, inscribedRadius, shape);
    else
        DUNE_THROW(Dune::NotImplemented, "Pc snap-off is not implemented for this irregular shape");
}

/*! \brief The entry capillary pressure of a pore throat.
 *
 * For details, see Eq. 11 in Rabbani et al., 2016: https://doi.org/10.1016/j.jcis.2016.03.053
 * or Eq A-7 in Oren et al., 1998: https://doi.org/10.2118/52052-PA
 */
template<class Scalar>
Scalar pcEntry(const Scalar surfaceTension,
               const Scalar contactAngle,
               const Scalar inscribedRadius,
               const Scalar shapeFactor) noexcept
{
    using std::sin;
    using std::cos;
    using std::sqrt;

    const Scalar cosContactAngle = cos(contactAngle);
    const Scalar sinContactAngle = sin(contactAngle);

    const Scalar D = M_PI - 3*contactAngle + 3*sinContactAngle*cosContactAngle - (cosContactAngle*cosContactAngle) / (4*shapeFactor);
    const Scalar F = (1 + sqrt(1 + 4*shapeFactor*D / (cosContactAngle*cosContactAngle))) / (1 + 2*sqrt(M_PI*shapeFactor));
    return surfaceTension / inscribedRadius * cosContactAngle * (1 + 2*sqrt(M_PI*shapeFactor)) * F;
}

} // namespace Dumux::PoreNetwork::ThresholdCapillaryPressures

#endif
