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
 * \brief Reation for a simple effective thermal conductivity
 */
#ifndef DUMUX_MATERIAL_THERMALCONDUCTIVITY_AVERAGE_HH
#define DUMUX_MATERIAL_THERMALCONDUCTIVITY_AVERAGE_HH

#include <algorithm>

namespace Dumux {

/*!
 * \ingroup Fluidmatrixinteractions
 * \brief Relation for a simple effective thermal conductivity
 */
template<class Scalar>
class ThermalConductivityAverage
{
public:
    /*!
     * \brief Relation for a simple effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
     *
     * \param volVars volume variables
     * \return Effective thermal conductivity \f$\mathrm{[W/(m K)]}\f$
     */
    template<class VolumeVariables>
    static Scalar effectiveThermalConductivity(const VolumeVariables& volVars)
    {
        constexpr int numFluidPhases = VolumeVariables::numFluidPhases();

        // Get the thermal conductivities and the porosity from the volume variables
        Scalar lambdaFluid = 0.0;
        for (int phaseIdx = 0; phaseIdx < numFluidPhases; ++phaseIdx)
            lambdaFluid += volVars.fluidThermalConductivity(phaseIdx)*volVars.saturation(phaseIdx);

        const Scalar lambdaSolid = volVars.solidThermalConductivity();
        const Scalar porosity = volVars.porosity();

        return lambdaSolid*(1-porosity) + lambdaFluid*porosity;
    }
};

} // end namespace Dumux

#endif
