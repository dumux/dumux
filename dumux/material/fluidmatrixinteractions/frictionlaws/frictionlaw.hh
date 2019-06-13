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
 * \copydoc Dumux::FrictionLaw
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_HH

#include <algorithm>
#include <cmath>

namespace Dumux {
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the abstract base class for friction laws.
 *
 * The LET mobility model is used to limit the friction for small water depths.
 */

template <typename Scalar, typename NumEqVector>
class FrictionLaw
{
public:
    /*!
     * \brief Compute the friction source term.
     *
     * The friction source term contains the losses due to friction.
     *
     * \return Friction source term.
     */

    virtual NumEqVector computeSource(const Scalar waterDept,
                                 const Scalar frictionValue,
                                 const Scalar u,
                                 const Scalar v) const = 0;

    /*!
     * \brief Limit the friction for small water depth.
     *
     * We define a water depth minUpperH. If the water depth is
     * smaller, we start to limit the friciton.
     * So the friciton term get's not extreme large for small water
     * depths.
     *
     * ------------------------- minUpperH -----------
     *
     *
     *
     * ------------------------roughnessHeight ---------------
     *    /\  /\   roughness                  /grain\
     * -------------------------------bottom ------------------
     * /////////////////////////////////////////////////
     *
     * For the limiting the LET model is used, which is usually applied in the
     * porous media flow to limit the permeability due to the saturation. It employs
     * the three empirical paramaters L, E and T, which describe the limiting curve (mobility).
     *
     * auto mobility = (mobility_max * pow(sw,L))/(pow(sw,L) + E * pow(1.0-sw,T));
     *
     * For the limitation of the roughness height L = 0.0, T = 2.0 and E = 1.0 are choosen.
     * Therefore the calculation of the mobility is simplified significantly.
     *
     * \param roughnessHeight roughness height of the representive structure (e.g. largest grain size).
     * \param waterDepth water depth.
     */
    Scalar limitRoughH(const Scalar roughnessHeight, const Scalar waterDepth) const
    {
        using std::min;
        using std::max;

        Scalar mobilityMax = 1.0; //!< maximal mobility

        Scalar minUpperH = roughnessHeight * 2.0;
        Scalar sw = min(waterDepth * (1.0/minUpperH),1.0);
        sw = max(0.0,sw);
        auto mobility = mobilityMax /(1 + (1.0-sw)*(1.0-sw));
        return roughnessHeight * (1.0 - mobility);
    }
};

} // end namespace Dumux

#endif // DUMUX_FRICTIONLAW_HH
