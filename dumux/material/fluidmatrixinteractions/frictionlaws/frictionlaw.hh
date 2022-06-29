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

#include <dune/common/fvector.hh>

namespace Dumux {
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the abstract base class for friction laws.
 *
 * The LET mobility model is used to limit the friction for small water depths.
 */

template <typename VolumeVariables >
class FrictionLaw
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:
    /*!
     * \brief Compute the bottom shear stress.
     *
     * \param volVars Volume variables
     *
     * Compute the bottom shear stress due to bottom friction.
     * The bottom shear stress is a projection of the shear stress tensor onto the river bed.
     * It can therefore be represented by a (tangent) vector with two entries.
     *
     * \return shear stress [N/m^2]. First entry is the x-component, the second the y-component.
     */
    virtual Dune::FieldVector<Scalar, 2> bottomShearStress(const VolumeVariables& volVars) const = 0;

    /*!
     * \brief Limit the friction for small water depth.
     *
     * We define a water depth minUpperH. If the water depth is
     * smaller, we start to limit the friction.
     * So the friction term get's not extreme large for small water
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
     * \param roughnessHeight roughness height of the representative structure (e.g. largest grain size).
     * \param waterDepth water depth.
     */
    Scalar limitRoughH(const Scalar roughnessHeight, const Scalar waterDepth) const
    {
        using std::clamp;

        const Scalar mobilityMax = 1.0; //!< maximal mobility

        const Scalar minUpperH = roughnessHeight * 2.0;
        const Scalar sw = clamp(waterDepth * (1.0/minUpperH), 0.0, 1.0);
        const Scalar mobility = mobilityMax /(1.0 + (1.0-sw)*(1.0-sw));
        return roughnessHeight * (1.0 - mobility);
    }

    // virtual base class needs a virtual destructor
    virtual ~FrictionLaw() = default;
};

} // end namespace Dumux

#endif
