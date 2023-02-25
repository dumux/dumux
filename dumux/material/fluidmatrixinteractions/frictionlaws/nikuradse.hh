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
#include <dune/common/math.hh>
#include "frictionlaw.hh"

namespace Dumux {
/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Implementation of the friction law after Nikuradse.
 *
 * The LET mobility model is used to limit the friction for small water
 * depths if a roughness height > 0.0 is provided (default roughnessHeight = 0.0).
 */

template <typename VolumeVariables>
class FrictionLawNikuradse : public FrictionLaw<VolumeVariables>
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:
    /*!
     * \brief Constructor
     *
     * \param ks Equivalent sand roughness (in m)
     * \param roughnessHeight roughness height (in m) default = 0.0
     */
    FrictionLawNikuradse(const Scalar ks, const Scalar roughnessHeight=0.0)
    : ks_(ks), roughnessHeight_(roughnessHeight) {}

    /*!
     * \brief Compute the bottom shear stress.
     *
     * \param volVars Volume variables
     *
     * Compute the bottom shear stress due to bottom friction.
     * The bottom shear stress is a projection of the shear stress tensor onto the river bed.
     * It can therefore be represented by a (tangent) vector with two entries.
     *
     * \return shear stress in N/m^2. First entry is the x-component, the second the y-component.
     */
    Dune::FieldVector<Scalar, 2> bottomShearStress(const VolumeVariables& volVars) const final
    {
        using Dune::power;
        using std::log;
        using std::hypot;

        Dune::FieldVector<Scalar, 2> shearStress(0.0);

        const Scalar artificialWaterDepth = this->limitRoughH(roughnessHeight_, volVars.waterDepth());
        const Scalar karmanConstant = 0.41; // Karman's constant is dimensionless
        const Scalar dimensionlessFactor = power(karmanConstant, 2)/power(log((12*(volVars.waterDepth() + artificialWaterDepth))/ks_), 2);
        const Scalar uv = hypot(volVars.velocity(0),volVars.velocity(1));

        shearStress[0] = dimensionlessFactor * volVars.velocity(0) * uv * volVars.density();
        shearStress[1] = dimensionlessFactor * volVars.velocity(1) * uv * volVars.density();

        return shearStress;
    }

private:
    Scalar ks_;
    Scalar roughnessHeight_;
};

} // end namespace Dumux

#endif
