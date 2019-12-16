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
 * \copydoc Dumux::FrictionLawNoFriction
 */
#ifndef DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_NOFRICTION_HH
#define DUMUX_MATERIAL_FLUIDMATRIX_FRICTIONLAW_NOFRICTION_HH

#include "frictionlaw.hh"

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief A pseudo friction law with no bottom friction
 */
template <typename VolumeVariables>
class FrictionLawNoFriction : public FrictionLaw<VolumeVariables>
{
    using Scalar = typename VolumeVariables::PrimaryVariables::value_type;
public:
    /*!
     * \brief Constructor
     */
    FrictionLawNoFriction() = default;

    /*!
     * \brief Compute the shear stress.
     *
     * \param volVars Volume variables
     *
     * Compute the shear stress due to friction. The shear stress is not a tensor as know
     * from contiuums mechanics, but a force projected on an area. Therefore it is a
     * vector with two entries. For this law without friction, the shearStress is zero.
     *
     * \return shear stress [N/m^2]. First entry is the x-component, the second the y-component.
     */
    Dune::FieldVector<Scalar, 2> shearStress(const VolumeVariables& volVars) const final
    {
        return {0.0, 0.0};
    }
};

} // end namespace Dumux

#endif
