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
 * \copydoc Dumux::FrictionLawManning
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_MANNING_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_MANNING_HH

#include <algorithm>
#include <cmath>
#include "frictionlaw.hh"

namespace Dumux {
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the friction law after Manning.
 *
 * The LET mobility model is used to limit the friction for small water depths.
 */

template <typename NumEqVector>
class FrictionLawManning : public FrictionLaw<NumEqVector>
{
    using Scalar = typename NumEqVector::value_type;

public:
    FrictionLawManning(const Scalar gravity)
        : gravity_(gravity) {}

    /*!
     * \brief Compute the friction source term.
     *
     * \param waterDepth water depth.
     * \param frictionValue Manning friction value.
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
        using std::hypot;

        NumEqVector source(0.0);

        Scalar roughnessHeight = pow(25.68/(1.0/frictionValue),6.0);
        roughnessHeight = this->limitRoughH(roughnessHeight, waterDepth);
        const Scalar c = pow((waterDepth + roughnessHeight),1.0/6.0) * 1.0/(frictionValue);
        const Scalar uv = hypot(u,v);

        source[0] = 0.0;
        source[1] = -gravity_/(c*c) * u * uv;
        source[2] = -gravity_/(c*c) * v * uv;

        return source;
    }
private:
    Scalar gravity_;
};

} // end namespace Dumux

#endif // DUMUX_FRICTIONLAW_MANNING_HH
