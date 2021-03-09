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
#include <dune/common/math.hh>
#include "frictionlaw.hh"

namespace Dumux {
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the friction law after Manning.
 *
 * The LET mobility model is used to limit the friction for small water depths.
 */

template <typename VolumeVariables>
class FrictionLawManning : public FrictionLaw<VolumeVariables>
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:
    /*!
     * \brief Constructor
     *
     * \param gravity Gravity constant [m/s^2]
     * \param manningN Manning friction coefficient [-]
     */
    FrictionLawManning(const Scalar gravity, const Scalar manningN)
        : gravity_(gravity), manningN_(manningN) {}

    /*!
     * \brief Compute the shear stress.
     *
     * \param volVars Volume variables
     *
     * Compute the shear stress due to friction. The shear stress is not a tensor as know
     * from contiuums mechanics, but a force projected on an area. Therefore it is a
     * vector with two entries.
     *
     * \return shear stress [N/m^2]. First entry is the x-component, the second the y-component.
     */
    Dune::FieldVector<Scalar, 2> shearStress(const VolumeVariables& volVars) const final
    {
        using std::pow;
        using Dune::power;
        using std::hypot;

        Dune::FieldVector<Scalar, 2> shearStress(0.0);

        Scalar roughnessHeight = power(25.68/(1.0/manningN_),6);
        roughnessHeight = this->limitRoughH(roughnessHeight, volVars.waterDepth());
        const Scalar c = pow((volVars.waterDepth() + roughnessHeight),1.0/6.0) * 1.0/(manningN_);
        const Scalar uv = hypot(volVars.velocity(0),volVars.velocity(1));

        shearStress[0] = -gravity_/(c*c) * volVars.velocity(0) * uv;
        shearStress[1] = -gravity_/(c*c) * volVars.velocity(1) * uv;

        return shearStress;
    }
private:
    Scalar gravity_;
    Scalar manningN_;
};

} // end namespace Dumux

#endif
