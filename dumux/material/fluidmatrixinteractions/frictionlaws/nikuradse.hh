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
 * \copydoc Dumux::FrictionLawNikuradse
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_NIKURADSE_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_NIKURADSE_HH

#include <algorithm>
#include <cmath>
#include "frictionlaw.hh"

namespace Dumux {
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the friction law after Nikuradse.
 *
 * The LET mobility model is used to limit the friction for small water depths.
 */

template <typename NumEqVector>
class FrictionLawNikuradse : public FrictionLaw<NumEqVector>
{
    using Scalar = typename NumEqVector::value_type;

public:
    /*!
     * \brief Compute the friction source term.
     *
     * \param waterDepth water depth.
     * \param frictionValue The equivalent sand roughness.
     * \param u velocity in x-direction.
     * \param v velocity in y-direction.
     *
     * \return Friction source term.
     */

    NumEqVector computeSource(const Scalar waterDepth,
                              const Scalar frictionValue,
                              const Scalar u,
                              const Scalar v) const final
    {
        using std::pow;
        using std::log;
        using std::hypot;

        NumEqVector source(0.0);

        Scalar roughnessHeight = frictionValue;
        roughnessHeight = this->limitRoughH(roughnessHeight, waterDepth);
        const Scalar ustarH = pow(0.41,2.0)/pow(log((12*(waterDepth + roughnessHeight))/frictionValue),2.0);
        const Scalar uv = hypot(u,v);

        source[0] = 0.0;
        source[1] = -ustarH * u * uv;
        source[2] = -ustarH * v * uv;

        return source;
    }
};

} // end namespace Dumux

#endif // DUMUX_FRICTIONLAW_NIKURADSE_HH
